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
  The screening kernel, the certification and the descriptive values all live in binom_core.py,
  which is the single implementation shared with dump_ties.py.  For each tie point (i,j) the masses
  come from the recurrence f_{k+1} = f_k * rho * (n-k)/(k+1) outward from the mode, the ranks from a
  two-pointer merge of the increasing left run and the decreasing right run (unimodality), then S_-
  and S_+.  No sorting, O(n) per tie point, O(n) memory.  Masses below TINY are set to zero.
  The masses are always normalised by their own sum; see binom_core.py for why.
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
from multiprocessing import Pool, cpu_count

import binom_core as core
from binom_core import MARGIN, GAP, TINY, TAG_MIN

HEADER = "n,i,j,pstar,E,F3,F3_sign,S_minus,S_plus,slope_left,slope_right,certified_by\n"

def screen_n(n):
    """[(i,j,'MIN'|'CHECK')] for every tie point of n that is not decided NOT in double precision."""
    r = core.screen(n, collect_all=False)
    return [(int(r['i'][t]), int(r['j'][t]), 'MIN' if r['tag'][t] == TAG_MIN else 'CHECK')
            for t in range(len(r['i']))]

certify = core.certify
evaluate = core.evaluate
recheck = core.recheck

def work(args):
    n, outdir = args
    path = os.path.join(outdir, f"n{n:05d}.csv")
    if os.path.exists(path): return n, -1, 0.0
    t0 = time.time(); rows = []; checks = []
    for (i, j, tag) in screen_n(n):
        how = 'double'
        if tag == 'CHECK':
            v, how = core.certify_escalating(n, i, j)
            if v is None: how = 'double'
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
    # strict name match: a sync service can leave stale copies like "n01993 2.csv" next to the
    # real files, and a loose n*.csv glob merges them in as silent duplicate rows.
    import re
    pat = re.compile(r"^n(\d{5})\.csv$")
    files = []
    for f in sorted(glob.glob(os.path.join(outdir, "n*.csv"))):
        m = pat.match(os.path.basename(f))
        if not m:
            print(f"  skipping unexpected file: {os.path.basename(f)}"); continue
        n = int(m.group(1))
        if (nmin is None or n >= nmin) and (nmax is None or n <= nmax): files.append(f)
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
