"""
cusps_fast.py -- cusp points (tie-point local minima) of the ordered-binomial expectation E(n,p),
for p* > 1/2 and n up to NMAX.  Parallel over n, resumable (one file per n).

    pip install numpy mpmath numba
    python cusps_fast.py --nmax 2000 --workers 8 --out cusps/
    python cusps_fast.py --merge --out cusps/                      # -> cusps/cusps_all.csv
    python cusps_fast.py --merge --out cusps/ --split-mb 25 --merge-nmin 1001
    python cusps_fast.py --recheck 34 12 30 --dps 50               # one row at high precision

Definitions
  f_p(k)  = C(n,k) p^k (1-p)^(n-k),  k = 0..n
  p*      = tie point of masses i<j (0<i<j<n, i+j>n so p*>1/2):  f(i)=f(j);
            rho = p*/(1-p*) = (C(n,i)/C(n,j))^(1/(j-i))
  w_k     = rank of f(k) in increasing order (0 = smallest).  Left of p*: w_j = w_i - 1.
  E(n,p)  = sum_k w_k f_p(k)
  S_-     = sum_k w_k^- f(k)(k - n p*)        (numerator of the left slope, E'_- = S_-/(p q))
  S_+     = S_- + (j-i) f(i)                  (numerator of the right slope, E'_+ = S_+/(p q))
  cusp    <=>  S_- < 0 < S_+   (local minimum of E at p*)
  F3      = (n+i-j)(i+j-2np*) + (j-np*)      (question (3) of the note)

Method
  Screening kernel (numba): for each tie point (i,j), masses via the recurrence
  f_{k+1} = f_k * rho * (n-k)/(k+1) outward from the mode, ranks via a two-pointer merge of the
  increasing left run and the decreasing right run (unimodality), then S_- and S_+.
  No sorting, O(n) per tie point, O(n) memory.  Masses below TINY are set to zero.
  A tie point is accepted in double precision if S_- < -MARGIN and S_+ > MARGIN (or rejected if
  clearly the other way).  It is instead re-certified with 50/100/200-digit interval arithmetic
  (mpmath.iv) if it is within MARGIN of a decision boundary, or if two adjacent masses in the
  ranking (both >= TINY) are within relative GAP of each other.  Every such tie point is logged in
  <out>/interval_checks.log with its verdict (MIN / NOT / UNRESOLVED) and the precision used.
  Descriptive columns (E, F3, slopes) are reported in double precision.

Output columns (one CSV per n, merged into cusps_all.csv)
  n,i,j,pstar,E,F3,F3_sign,S_minus,S_plus,slope_left,slope_right,certified_by
  certified_by = double | iv50 | iv100 | iv200;  an undecided tie point is written as UNRESOLVED.

History
  This replaces an earlier numpy-only script (cusps_parallel.py) that ranked the masses with a
  full sort and flagged a near-tie whenever ANY two adjacent sorted masses were within relative
  1e-8, including deep-tail masses (~1e-300) that cannot affect the sign of S_-/S_+.  It was fine
  up to n~1080 (12 interval checks at n=1000) but the count exploded beyond that (182,542 checks,
  ~12 h, at n=2000).  Ignoring masses below TINY brings n=2000 down to ~75 checks and ~40 s.
  Verified: rerunning n=3..1000 with this script reproduces every per-n file and the
  interval-check log of the old script byte for byte.  See README_cusps.md.
"""
import argparse, os, sys, time, glob
import numpy as np
from math import comb, lgamma
from multiprocessing import Pool, cpu_count
from numba import njit

MARGIN, GAP, TINY = 1e-6, 1e-8, 1e-290
HEADER = "n,i,j,pstar,E,F3,F3_sign,S_minus,S_plus,slope_left,slope_right,certified_by\n"

@njit(cache=True)
def _screen_kernel(n, lnC, out_i, out_j, out_tag):
    """Fill out_* with tie points that are MIN (tag 1) or need CHECK (tag 2). Returns count."""
    cap = out_i.shape[0]
    f = np.empty(n+1); w = np.empty(n+1, np.int64); cnt = 0
    for i in range(1, n):
        for j in range(max(i+1, n-i+1), n):          # i<j<n and i+j>n
            m = j - i
            lnrho = (lnC[i] - lnC[j]) / m
            p = 1.0/(1.0 + np.exp(-lnrho)); q = 1.0 - p; rho = p/q
            # mode of Bin(n,p): floor((n+1)p)
            md = int(np.floor((n+1)*p))
            if md > n: md = n
            f[md] = np.exp(lnC[md] + md*np.log(p) + (n-md)*np.log(q))
            for k in range(md, 0, -1):               # leftwards: f_{k-1} = f_k * k/((n-k+1) rho)
                f[k-1] = f[k] * k / ((n-k+1.0)*rho)
            for k in range(md, n):                   # rightwards
                f[k+1] = f[k] * rho*(n-k) / (k+1.0)
            for k in range(n+1):                     # masses below TINY are numerically zero
                if f[k] < TINY: f[k] = 0.0
            f[j] = f[i]                              # exact tie
            lnkap = np.log(m) + lnC[i] + i*np.log(p) + (n-i)*np.log(q)   # kink height (j-i) f(i), exact
            kap = np.exp(lnkap) if lnkap > -700.0 else 0.0
            # split point: left run 0..s increasing, right run s+1..n decreasing; need i<=s<j
            s = md
            if s >= j: s = j - 1
            if s < i: s = i
            # two-pointer merge, increasing order; ties: j before i, otherwise left first
            a = 0; b = n; r = 0
            neartie = False; prev = -1.0
            while a <= s or b > s:
                if a > s: take_left = False
                elif b <= s: take_left = True
                elif f[a] < f[b]: take_left = True
                elif f[a] > f[b]: take_left = False
                else:                                # equal
                    take_left = not (a == i and b == j)
                if take_left: k = a; a += 1
                else:         k = b; b -= 1
                w[k] = r
                if r > 0 and f[k] > 0.0 and not ((k == i and prev == f[j]) or (k == j and prev == f[i])):
                    if (f[k] - prev)/f[k] < GAP and not (k == i or k == j): neartie = True
                    if (k == i or k == j) and (f[k]-prev)/f[k] < GAP and prev != f[k]: neartie = True
                prev = f[k]; r += 1
            Sm = 0.0
            for k in range(n+1): Sm += w[k]*f[k]*(k - n*p)
            Sp = Sm + kap
            ismin = (Sm < -MARGIN) and (Sp > MARGIN)
            isnot = (Sm > MARGIN) or (Sp < -MARGIN)
            if neartie or not (ismin or isnot):
                if cnt < cap:
                    out_i[cnt] = i; out_j[cnt] = j; out_tag[cnt] = 2
                cnt += 1
            elif ismin:
                if cnt < cap:
                    out_i[cnt] = i; out_j[cnt] = j; out_tag[cnt] = 1
                cnt += 1
    return cnt

def screen_n(n):
    lnC = np.array([lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1) for k in range(n+1)])
    cap = max(64, 4*n)
    while True:
        oi = np.empty(cap, np.int64); oj = np.empty(cap, np.int64); ot = np.empty(cap, np.int64)
        c = _screen_kernel(n, lnC, oi, oj, ot)
        if c <= cap: break
        cap = 2*c                                    # too small: retry with a bigger buffer
    return [(int(oi[t]), int(oj[t]), 'MIN' if ot[t] == 1 else 'CHECK') for t in range(c)]

# --- certification, descriptive values, per-n driver, merge, recheck -------------------
def certify(n, i, j, dps):
    from mpmath import iv
    iv.dps = dps; m = j-i
    rho = (iv.mpf(comb(n,i))/iv.mpf(comb(n,j)))**(iv.mpf(1)/m)
    p = rho/(1+rho); q = 1-p
    f = [iv.mpf(comb(n,k))*p**k*q**(n-k) for k in range(n+1)]
    order = sorted(range(n+1), key=lambda k: (f[i].mid if k in (i,j) else f[k].mid, 0 if k==j else 1))
    for a, b in zip(order, order[1:]):
        if {a,b} == {i,j}: continue
        if not (f[a].b < f[b].a): return None
    w = [0]*(n+1)
    for r, k in enumerate(order): w[k] = r
    Sm = sum(w[k]*f[k]*(k-n*p) for k in range(n+1)); Sp = Sm + m*f[i]
    if Sm.b < 0 and Sp.a > 0: return 'MIN'
    if Sm.a >= 0 or Sp.b <= 0: return 'NOT'
    return None

def evaluate(n, i, j):
    lnC = np.array([lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1) for k in range(n+1)]); k = np.arange(n+1); m = j-i
    lnrho = (lnC[i]-lnC[j])/m; p = 1/(1+np.exp(-lnrho)); q = 1-p
    with np.errstate(under='ignore'):
        f = np.exp(lnC + k*np.log(p) + (n-k)*np.log(q))
    f[j] = f[i]
    key2 = np.zeros(n+1); key2[i] = 1
    order = np.lexsort((key2, f)); w = np.empty(n+1, int); w[order] = k
    Sm = float(np.sum(w*f*(k-n*p))); Sp = Sm + m*f[i]
    E = float(np.sum(w*f)); F3 = (n+i-j)*(i+j-2*n*p) + (j-n*p)
    return p, E, F3, Sm, Sp, Sm/(p*q), Sp/(p*q)

def work(args):
    n, outdir = args
    path = os.path.join(outdir, f"n{n:05d}.csv")
    if os.path.exists(path): return n, -1, 0.0
    t0 = time.time(); rows = []; checks = []
    for (i, j, tag) in screen_n(n):
        how = 'double'
        if tag == 'CHECK':
            v = None
            for dps in (50, 100, 200):
                v = certify(n, i, j, dps)
                if v: how = f'iv{dps}'; break
            checks.append(f"{n},{i},{j},{v or 'UNRESOLVED'},{how}\n")
            if v != 'MIN':
                if v is None: rows.append((i, j, 'UNRESOLVED'))
                continue
        rows.append((i, j, how))
    if checks:
        with open(os.path.join(outdir, "interval_checks.log"), 'a') as lg: lg.writelines(checks)
    with open(path + '.tmp', 'w') as fo:
        fo.write(HEADER)
        for (i, j, how) in rows:
            if how == 'UNRESOLVED':
                fo.write(f"{n},{i},{j},,,,,,,,,UNRESOLVED\n"); continue
            p, E, F3, Sm, Sp, sl, sr = evaluate(n, i, j)
            fo.write(f"{n},{i},{j},{p:.16g},{E:.16g},{F3:.16g},{'+' if F3>0 else '-'},"
                     f"{Sm:.16g},{Sp:.16g},{sl:.16g},{sr:.16g},{how}\n")
    os.replace(path + '.tmp', path)
    return n, len(rows), time.time()-t0

def merge(outdir, split_mb=0, nmin=None, nmax=None, slim=False):
    files = sorted(glob.glob(os.path.join(outdir, "n*.csv")))
    files = [f for f in files if (nmin is None or int(os.path.basename(f)[1:6]) >= nmin)
                              and (nmax is None or int(os.path.basename(f)[1:6]) <= nmax)]
    header = HEADER
    def conv(line):
        if not slim: return line
        c = line.rstrip("\n").split(",")
        if c[3] == "": return ",".join(c[:3]) + ",,,,,,\n"
        g = lambda x: f"{float(x):.12g}"
        return f"{c[0]},{c[1]},{c[2]},{g(c[3])},{g(c[4])},{g(c[5])},{c[6]},{g(c[9])},{g(c[10])}\n"
    if slim: header = "n,i,j,pstar,E,F3,F3_sign,slope_left,slope_right\n"
    limit = split_mb*1024*1024 if split_mb else None
    part = 0; fo = None; size = 0; manifest = []; total = 0; first_n = None; last_n = None
    def open_part():
        nonlocal part, fo, size, first_n
        part += 1; name = os.path.join(outdir, f"cusps_part{part:02d}.csv" if limit else "cusps_all.csv")
        fo = open(name, "w"); fo.write(header); size = len(header); first_n = None; return name
    name = open_part()
    for fn in files:
        n = int(os.path.basename(fn)[1:6])
        with open(fn) as fi:
            next(fi); lines = [conv(l) for l in fi]
        nbytes = sum(len(l) for l in lines)
        if limit and size + nbytes > limit and first_n is not None:
            fo.close(); manifest.append((name, first_n, last_n, size)); name = open_part()
        if first_n is None: first_n = n
        fo.writelines(lines); size += nbytes; total += len(lines); last_n = n
    fo.close(); manifest.append((name, first_n, last_n, size))
    with open(os.path.join(outdir, "cusps_manifest.txt"), "w") as mf:
        for name, a, b, sz in manifest:
            line = f"{os.path.basename(name)}: n={a}..{b}, {sz/1e6:.1f} MB"; print(line); mf.write(line + "\n")
    print(f"{total} rows in {len(manifest)} file(s)")

def recheck(n, i, j, dps=50):
    from mpmath import mp, mpf, nstr
    mp.dps = dps; m = j-i
    rho = (mpf(comb(n,i))/comb(n,j))**(mpf(1)/m); p = rho/(1+rho); q = 1-p
    f = [comb(n,k)*p**k*q**(n-k) for k in range(n+1)]
    order = sorted(range(n+1), key=lambda k: (f[i] if k in (i,j) else f[k], 0 if k==j else 1))
    w = [0]*(n+1)
    for r, k in enumerate(order): w[k] = r
    Sm = sum(w[k]*f[k]*(k-n*p) for k in range(n+1)); Sp = Sm + m*f[i]
    E = sum(w[k]*f[k] for k in range(n+1)); F3 = (n+i-j)*(i+j-2*n*p) + (j-n*p)
    print(f"n={n} i={i} j={j}  ({dps} digits)")
    for name, v in (("p*",p),("E",E),("F3",F3),("S_-",Sm),("S_+",Sp),("slope_left",Sm/(p*q)),("slope_right",Sp/(p*q))):
        print(f"  {name:12s} {nstr(v, dps-10)}")
    print("  cusp:", Sm < 0 < Sp)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=200)
    ap.add_argument("--nmin", type=int, default=3)
    ap.add_argument("--workers", default="auto")
    ap.add_argument("--out", default="cusps")
    ap.add_argument("--merge", action="store_true")
    ap.add_argument("--split-mb", type=float, default=0)
    ap.add_argument("--slim", action="store_true")
    ap.add_argument("--merge-nmin", type=int); ap.add_argument("--merge-nmax", type=int)
    ap.add_argument("--recheck", nargs=3, type=int, metavar=("N","I","J"))
    ap.add_argument("--dps", type=int, default=50)
    a = ap.parse_args()
    if a.recheck: recheck(*a.recheck, dps=a.dps); sys.exit()
    if a.merge: merge(a.out, a.split_mb, a.merge_nmin, a.merge_nmax, a.slim); sys.exit()
    os.makedirs(a.out, exist_ok=True)
    screen_n(10)                                     # compile the kernel once in the parent
    workers = cpu_count() if a.workers == "auto" else int(a.workers)
    ns = list(range(a.nmax, a.nmin-1, -1))
    weight = {n: n**3 for n in ns}; total_w = sum(weight.values())
    already = [n for n in ns if os.path.exists(os.path.join(a.out, f"n{n:05d}.csv"))]
    if already: print(f"resuming: {len(already)} of {len(ns)} values of n already done", flush=True)
    t0 = time.time(); done = 0; tot = 0; done_w = sum(weight[n] for n in already); last = t0
    print(f"screening n={a.nmin}..{a.nmax} with {workers} workers (largest n first)", flush=True)
    with Pool(workers) as pool:
        for n, cnt, dt in pool.imap_unordered(work, [(n, a.out) for n in ns]):
            done += 1
            if cnt >= 0: tot += cnt; done_w += weight[n]
            now = time.time()
            if now - last > 15 or done == len(ns):
                frac = done_w/total_w; el = now - t0
                eta = el*(1-frac)/max(frac, 1e-9) if frac < 1 else 0
                print(f"[{el:6.0f}s] {done:4d}/{len(ns)} n-values | {100*frac:5.1f}% of work | "
                      f"last n={n} ({cnt} cusps, {dt:.1f}s) | cusps so far {tot} | ETA {eta/60:.1f} min", flush=True)
                last = now
    print(f"finished in {time.time()-t0:.0f}s with {workers} workers. Now run:  python {sys.argv[0]} --merge --out {a.out}")
