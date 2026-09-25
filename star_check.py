"""
star_check.py -- exhaustive switch-point check of the (★) reduction, per n, in parallel.

    .venv/bin/python star_check.py --validate                    # table of RESEARCH_LOG 2026-09-24 item 4b
    .venv/bin/python star_check.py --mp-check                    # 50-digit spot checks (mpmath)
    .venv/bin/python star_check.py --lemma-check 300 --workers 8 # sampled check of the switch-point lemma
    .venv/bin/python star_check.py --nmax 5000 --workers 8 --out analysis/

Definitions (RESEARCH_LOG.md section 4 and the 2026-09-24 claude.ai entry, items 1-3)
  D(p)   = min over c in Z/2 + 1/4 of E_p|K - c|,  K ~ Bin(n,p).
  L(p)   = n + 1/2 - 2 D(p) <= E(n,p) for all p, with equality at p = 1/2, so
           (★)  D(p) <= D(1/2) for all p   ==>   E(n,p) >= E(n,1/2).
  Switch point p_m, m = 1..n-1, the m-th point above 1/2 where the minimising centre moves;
  a = floor((n+m)/2):
     HALF type (n+m odd):  P(K <= a) = 1/2;             centre a+1/4 -> a+3/4
     INT  type (n+m even): P(K < a) = P(K > a);         centre a-1/4 -> a+1/4
  (m = n would be INT with a = n, which has no root below p = 1, so there are exactly n-1.)
  By the lemma (item 2, ARGUED) D is quasi-convex between switch points, so (★) <=> D(p_m) <= D(1/2)
  for every m.  Scaled margin M(n,m) = 2 (D(1/2) - D(p_m)) n^1.5;  x_m = n (p_m - 1/2).

Method
  Roots: both balance functions are decreasing in p, positive at p = 1/2 (m >= 1) and negative at
  p = 1.  Each root is bracketed near 1/2 + m/(2n) (verified; fallback [1/2, 1]) and bisected,
  vectorised over m, until the bracket is two ADJACENT doubles; p_m is the end with the smaller
  |residual|.  That is the best any method can do with the computed residual -- no p-grid anywhere.
  Cdfs are scipy.special.bdtr/bdtrc (regularised incomplete beta).
  D without cancellation: with K' ~ Bin(n-1,p), F' its cdf and f' its pmf,
        E_p|K - a| = (a - np)(2 F'(a) - 1) + 2 p (n - a) f'(a)                 (exact identity)
  (from E|K-a| = np - a + 2[a P(K<=a) - np P(K'<=a-1)] and P(K<=a) = F'(a) - p f'(a)).
  Near the median a - np = O(1) and 2F'-1 = O(n^-1/2), so nothing of size n is subtracted.  Then
  the lower envelope over the two centres that meet at the switch, evaluated at the computed p:
        HALF:  D = E|K-a| + min(s/4, 3s/4),        s = 2 P(K<=a) - 1
        INT:   D = E|K-a| + (P(K=a) - |t|)/4,     t = P(K<a) - P(K>a)
  (s, t = 0 at the exact root, where these reduce to item 3's E|K-a-1/2| and E|K-a| + P(K=a)/4.)
  D(1/2) is EXACT: 4 * 2^n * E|K - (n/2 + 1/4)| is an integer, computed with Python integers and
  rounded once.
  Cross-check, every run: for a seeded random sample of (n,m), D(p_m) again by direct summation over
  k, minimised over EVERY centre c in Z/2 + 1/4 (prefix sums), from normalised lgamma masses.

Outputs (LF line endings)
  <out>/star_per_n.csv            n,n_switch,min_M,argmin_m,argmin_type,M1,M2,M3,M4,x1,D_half
  <out>/star_switch_points/nNNNNN.csv    n,m,type,a,p_m,x_m,D_pm,M   (local only: ~n^2 in total)
  <out>/star_crosscheck.csv       the direct-summation sample
  <out>/star_timings.csv          n,seconds  (wall time per n inside the worker)
"""
import argparse, csv, math, os, random, sys, time
from fractions import Fraction
from math import comb, lgamma
from multiprocessing import Pool, cpu_count

import numpy as np
from scipy.special import bdtr, bdtrc
from scipy.stats import binom

HALF, INT = "HALF", "INT"

# ---------------------------------------------------------------------------------- switch points
def balance(n, m, a, p):
    """Balance residual (decreasing in p; zero at p_m).  m, a, p are arrays."""
    half = (n + m) % 2 == 1
    g = np.empty_like(p)
    h = half
    g[h] = bdtr(a[h], n, p[h]) - 0.5
    g[~h] = bdtr(a[~h] - 1, n, p[~h]) - bdtrc(a[~h], n, p[~h])
    return g

def switch_points(n):
    """(m, a, is_half, p_m, residual) for m = 1..n-1, roots bisected to adjacent doubles."""
    m = np.arange(1, n)
    a = (n + m) // 2
    guess = 0.5 + m/(2.0*n)
    lo = np.maximum(0.5, guess - 1.0/n)
    hi = np.minimum(1.0, guess + 1.0/n)
    bad = ~((balance(n, m, a, lo) > 0) & (balance(n, m, a, hi) < 0))
    lo[bad] = 0.5; hi[bad] = 1.0
    while True:
        mid = 0.5*(lo + hi)
        act = (mid > lo) & (mid < hi)                   # False once lo, hi are adjacent doubles
        if not act.any(): break
        ga = balance(n, m[act], a[act], mid[act])
        idx = np.flatnonzero(act)
        pos = ga > 0
        lo[idx[pos]] = mid[act][pos]
        hi[idx[~pos]] = mid[act][~pos]
    glo = balance(n, m, a, lo); ghi = balance(n, m, a, hi)
    p = np.where(np.abs(glo) <= np.abs(ghi), lo, hi)
    res = np.where(np.abs(glo) <= np.abs(ghi), glo, ghi)
    return m, a, (n + m) % 2 == 1, p, res

def mad_int(n, a, p):
    """E_p|K - a| for integer a, cancellation-free (see docstring)."""
    F1 = bdtr(a, n - 1, p)
    f1 = binom.pmf(a, n - 1, p)
    return (a - n*p)*(2.0*F1 - 1.0) + 2.0*p*(n - a)*f1

def D_at(n, a, is_half, p):
    """Lower envelope of the two centres that meet at the switch, at the computed p."""
    base = mad_int(n, a, p)
    s = 2.0*bdtr(a, n, p) - 1.0
    t = bdtr(a - 1, n, p) - bdtrc(a, n, p)
    pk = binom.pmf(a, n, p)
    return np.where(is_half, base + np.minimum(0.25*s, 0.75*s), base + 0.25*(pk - np.abs(t)))

def D_half_exact(n):
    """D(1/2), exactly rational, rounded once.  Section 4: E_{1/2}|K - n/2| = (n/2) C(n-1, floor(n/2)) 2^(1-n),
    plus P(K = n/2)/4 for n even (INT type at 1/2); for n odd the centre n/2 +- 1/4 sees a flat
    E|K-c| on [(n-1)/2, (n+1)/2], so D(1/2) = E|K - n/2|."""
    D = Fraction(n*comb(n - 1, n//2), 2**n)
    if n % 2 == 0: D += Fraction(comb(n, n//2), 4 * 2**n)
    return float(D)

def D_half_bruteforce(n):
    """The same, as the definition: E_{1/2}|K - (n/2 + 1/4)| summed exactly (O(n^2) bits; checks only)."""
    num = sum(comb(n, k)*abs(4*k - 2*n - 1) for k in range(n + 1))     # 4c = 2n+1
    return float(Fraction(num, 4 * 2**n))

def D_direct(n, p):
    """min over every c in Z/2 + 1/4 of E_p|K-c|, by direct summation (normalised lgamma masses)."""
    k = np.arange(n + 1)
    lnC = np.array([lgamma(n+1) - lgamma(j+1) - lgamma(n-j+1) for j in range(n + 1)])
    f = np.exp(lnC + k*math.log(p) + (n - k)*math.log1p(-p))
    f = f / math.fsum(f.tolist())
    F = np.cumsum(f); G = np.cumsum(f*k)                  # sums over k <= j
    tot, totk = F[-1], G[-1]
    best = np.inf
    for off in (0.25, 0.75):
        c = k + off                                       # c = j + off, floor(c) = j
        below = c*F - G                                   # sum_{k<=j} (c-k) f
        above = (totk - G) - c*(tot - F)                  # sum_{k>j} (k-c) f
        best = min(best, float(np.min(below + above)))
    best = min(best, float(np.sum(f*(k + 0.25))))         # c = -1/4
    return best

# --------------------------------------------------------------------------------- one n (worker)
def run_n(n, out=None):
    t0 = time.time()
    m, a, is_half, p, res = switch_points(n)
    Dp = D_at(n, a, is_half, p)
    Dh = D_half_exact(n)
    M = 2.0*(Dh - Dp)*n**1.5
    x = n*(p - 0.5)
    mono = bool(np.all(np.diff(p) > 0)) and p[0] > 0.5 and p[-1] < 1.0
    typ = np.where(is_half, HALF, INT)
    if out is not None:
        path = os.path.join(out, "star_switch_points", f"n{n:05d}.csv")
        with open(path + ".tmp", "w", newline="") as fh:
            w = csv.writer(fh, lineterminator="\n")
            w.writerow(["n", "m", "type", "a", "p_m", "x_m", "D_pm", "M"])
            for r in zip(m.tolist(), typ.tolist(), a.tolist(), p.tolist(), x.tolist(),
                         Dp.tolist(), M.tolist()):
                w.writerow((n,) + r)
        os.replace(path + ".tmp", path)
    i = int(np.argmin(M))
    Ms = [repr(float(M[j])) if j < len(M) else "" for j in range(4)]
    row = [n, n - 1, repr(float(M[i])), int(m[i]), typ[i]] + Ms + [repr(float(x[0])), repr(Dh)]
    return row, mono, float(np.max(np.abs(res))), M[:20].copy(), time.time() - t0

def crosscheck(n, mm):
    """Direct-summation D at p_m for one (n, m); returns (D formula, D direct, rel diff)."""
    m, a, is_half, p, _ = switch_points(n)
    j = mm - 1
    Df = float(D_at(n, a[j:j+1], is_half[j:j+1], p[j:j+1])[0])
    Dd = D_direct(n, float(p[j]))
    return Df, Dd, abs(Df - Dd)/Dd

def work(args):
    n, out, checks = args
    row, mono, maxres, M20, dt = run_n(n, out)
    cc = [(n, mm) + crosscheck(n, mm) for mm in checks]
    return n, row, mono, maxres, M20, dt, cc

# ------------------------------------------------------------------------------------ validation
TABLE = {  # RESEARCH_LOG 2026-09-24 (claude.ai) item 4b: x -> M at n = 1000, 1001, 3000, 3001
    1000: {1: 0.2990, 2: 1.5948, 3: 3.4885, 4: 6.3790, 5: 9.8676, 6: 14.3529, 8: 25.5164},
    1001: {1: 0.4983, 2: 1.5948, 3: 3.6879, 4: 6.3790, 5: 10.0669, 6: 14.3529, 8: 25.5164},
    3000: {1: 0.2991, 2: 1.5954, 3: 3.4900, 4: 6.3817, 5: 9.8717, 6: 14.3589, 8: 25.5269},
    3001: {1: 0.4986, 2: 1.5954, 3: 3.6894, 4: 6.3817, 5: 10.0712, 6: 14.3589, 8: 25.5269},
}

def validate():
    from scipy.optimize import brentq
    ok = True
    for n, tab in TABLE.items():
        row, mono, maxres, M20, dt = run_n(n)
        for mm, want in tab.items():
            got = float(M20[mm - 1]); good = round(got, 4) == want
            ok &= good
            print(f"n={n:5d} m={mm}  M={got:.6f}  table {want:.4f}  {'ok' if good else 'MISMATCH'}")
        # brentq on the same balance function, as an independent root finder
        m, a, is_half, p, _ = switch_points(n)
        for mm in (1, 2, 3, 4, n // 2, n - 1):
            j = mm - 1
            g = lambda q: float(balance(n, m[j:j+1], a[j:j+1], np.array([q]))[0])
            pb = brentq(g, 0.5, 1.0, xtol=1e-18, rtol=8.9e-16, maxiter=500)
            print(f"   brentq m={mm:5d}: |p_bisect - p_brentq| = {abs(pb - p[j]):.2e}")
        print(f"   monotone={mono}  max|residual|={maxres:.2e}  time={dt:.3f}s")
    bad = [n for n in list(range(3, 301)) + [999, 1000, 1001, 2000, 4999, 5000]
           if D_half_exact(n) != D_half_bruteforce(n)]
    print(f"D(1/2) closed form == exact definitional sum, bit for bit, n=3..300 and 6 larger n: "
          f"{'yes' if not bad else bad}")
    ok &= not bad
    print("VALIDATION", "PASSED" if ok else "FAILED")
    return ok

# ------------------------------------------------------------------------------- mpmath spot check
MP_PAIRS = [(3, 1), (3, 2), (4, 3), (10, 5), (57, 1), (100, 1), (100, 50), (100, 99), (777, 3),
            (1000, 1), (1000, 4), (1001, 1), (1234, 20), (2000, 1), (2000, 1000), (2500, 17),
            (3001, 1), (3001, 2), (4000, 3999), (4999, 1), (5000, 1), (5000, 2), (5000, 3),
            (5000, 2500), (5000, 4999)]

def mp_point(n, mm, p0, dps=50):
    """50-digit p_m, D(p_m), D(1/2) by direct summation.  Returns mpf values."""
    from mpmath import mp, mpf, findroot
    mp.dps = dps
    a = (n + mm)//2; half = (n + mm) % 2 == 1
    def masses(p):
        q = 1 - p; r = p/q
        f = [q**n]
        for k in range(n): f.append(f[-1]*(n - k)/(k + 1)*r)
        return f
    def g(p):
        f = masses(p)
        if half: return mp.fsum(f[:a + 1]) - mpf(1)/2
        return mp.fsum(f[:a]) - mp.fsum(f[a + 1:])
    eps = mpf(10)**-12
    pm = findroot(g, (mpf(p0) - eps, mpf(p0) + eps), solver="anderson", tol=mpf(10)**(-2*dps + 10))
    f = masses(pm)
    def mad(c): return mp.fsum(fk*abs(k - c) for k, fk in enumerate(f))
    cs = [a + mpf(o)/4 for o in (-3, -1, 1, 3, 5)]
    vals = [mad(c) for c in cs]
    j = min(range(5), key=lambda t: vals[t])
    assert 0 < j < 4, (n, mm, j)
    Dh = mp.fsum(comb(n, k)*abs(k - mpf(n)/2 - mpf(1)/4) for k in range(n + 1)) / mpf(2)**n
    return pm, vals[j], Dh

def mp_check(out=None):
    from mpmath import mp, mpf
    rows = []
    print(f"{'n':>5} {'m':>5} {'type':>4} {'|dp|/ulp':>9} {'|dD_pm|':>9} {'|dD_half|':>9} "
          f"{'|dD| at mp p':>12} {'M_mp':>12} {'|dM|':>9} {'rel dM':>9}")
    for n, mm in MP_PAIRS:
        m, a, is_half, p, _ = switch_points(n)
        j = mm - 1
        pd = float(p[j]); Dd = float(D_at(n, a[j:j+1], is_half[j:j+1], p[j:j+1])[0])
        Dh = D_half_exact(n); Md = 2*(Dh - Dd)*n**1.5
        pm, Dm, Dhm = mp_point(n, mm, pd)
        Mm = 2*(Dhm - Dm)*mpf(n)**mpf(1.5)
        pr = float(pm)                                           # mp root rounded to double
        Dr = float(D_at(n, a[j:j+1], is_half[j:j+1], np.array([pr]))[0])
        r = dict(n=n, m=mm, type=HALF if is_half[j] else INT,
                 dp_ulp=float(abs(pd - pm))/math.ulp(pr), dD_pm=float(abs(Dd - Dm)),
                 dD_half=float(abs(Dh - Dhm)), dD_at_mp_p=float(abs(Dr - Dm)),
                 D_pm=float(Dm), M_mp=float(Mm), dM=float(abs(Md - Mm)), rel_dM=float(abs(Md - Mm)/abs(Mm)))
        rows.append(r)
        print(f"{n:5d} {mm:5d} {r['type']:>4} {r['dp_ulp']:9.2f} {r['dD_pm']:9.2e} {r['dD_half']:9.2e} "
              f"{r['dD_at_mp_p']:12.2e} {r['M_mp']:12.6f} {r['dM']:9.2e} {r['rel_dM']:9.2e}", flush=True)
    if out:
        with open(os.path.join(out, "star_mp_check.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
            w.writeheader(); w.writerows(rows)
    return rows

# ------------------------------------------------------------------------------------ lemma check
def lemma_n(n, per=64):
    """Sample D (direct, all centres) inside every switch interval; report the worst excess of D over
    the larger endpoint value, and over D(1/2).  Quasi-convexity predicts excess <= 0."""
    m, a, is_half, p, _ = switch_points(n)
    Dp = D_at(n, a, is_half, p)
    ends = np.concatenate([[0.5], p, [1.0]])
    Dend = np.concatenate([[D_half_exact(n)], Dp, [0.25]])
    worst = -np.inf; worst_half = -np.inf
    for lo, hi, dlo, dhi in zip(ends[:-1], ends[1:], Dend[:-1], Dend[1:]):
        for t in (np.arange(1, per + 1) / (per + 1)):
            q = lo + t*(hi - lo)
            d = D_direct(n, q)
            worst = max(worst, d - max(dlo, dhi))
            worst_half = max(worst_half, d - Dend[0])
    return n, worst, worst_half

# -------------------------------------------------------------------------------------------- main
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=3)
    ap.add_argument("--nmax", type=int, default=5000)
    ap.add_argument("--workers", default="auto")
    ap.add_argument("--out", default="analysis")
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--mp-check", action="store_true")
    ap.add_argument("--lemma-check", type=int, metavar="NMAX")
    ap.add_argument("--n-cross", type=int, default=200, help="random (n,m) direct-summation checks")
    ap.add_argument("--seed", type=int, default=20260924)
    a = ap.parse_args()
    workers = cpu_count() if a.workers == "auto" else int(a.workers)
    if a.validate: sys.exit(0 if validate() else 1)
    if a.mp_check: mp_check(a.out); sys.exit()
    if a.lemma_check:
        t0 = time.time(); W = []; WH = []
        with Pool(workers) as pool:
            for n, w_, wh in pool.imap_unordered(lemma_n, range(a.lemma_check, 2, -1)):
                W.append((w_, n)); WH.append((wh, n))
        print(f"lemma check n=3..{a.lemma_check}, 64 samples per switch interval, {time.time()-t0:.0f}s")
        print(f"  max over n of [D(sample) - max(D at the interval's ends)] = {max(W)[0]:.3e} (n={max(W)[1]})")
        print(f"  max over n of [D(sample) - D(1/2)]                        = {max(WH)[0]:.3e} (n={max(WH)[1]})")
        sys.exit()

    os.makedirs(os.path.join(a.out, "star_switch_points"), exist_ok=True)
    ns = list(range(a.nmax, a.nmin - 1, -1))
    rng = random.Random(a.seed)
    checks = {n: [] for n in ns}
    for _ in range(a.n_cross):
        n = rng.randint(a.nmin, a.nmax); checks[n].append(rng.randint(1, n - 1))
    for n in (a.nmax, a.nmin):
        checks[n].append(1)
    t0 = time.time(); rows = {}; times = {}; cross = []; nonmono = []; maxres = 0.0
    total_w = sum(ns); done_w = 0; last = t0
    print(f"(★) switch-point check n={a.nmin}..{a.nmax} on {workers} workers", flush=True)
    with Pool(workers) as pool:
        for n, row, mono, mr, M20, dt, cc in pool.imap_unordered(
                work, [(n, a.out, checks[n]) for n in ns]):
            rows[n] = row; times[n] = dt; cross += cc; maxres = max(maxres, mr)
            if not mono: nonmono.append(n)
            done_w += n; now = time.time()
            if now - last > 15:
                print(f"[{now-t0:6.0f}s] {100*done_w/total_w:5.1f}% of work, last n={n} ({dt:.2f}s)", flush=True)
                last = now
    wall = time.time() - t0
    with open(os.path.join(a.out, "star_per_n.csv"), "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["n", "n_switch", "min_M", "argmin_m", "argmin_type", "M1", "M2", "M3", "M4", "x1", "D_half"])
        for n in sorted(rows): w.writerow(rows[n])
    with open(os.path.join(a.out, "star_timings.csv"), "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n"); w.writerow(["n", "seconds"])
        for n in sorted(times): w.writerow([n, f"{times[n]:.4f}"])
    with open(os.path.join(a.out, "star_crosscheck.csv"), "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n"); w.writerow(["n", "m", "D_formula", "D_direct", "rel_diff"])
        for r in sorted(cross): w.writerow([r[0], r[1], repr(r[2]), repr(r[3]), f"{r[4]:.3e}"])
    minM = min((float(r[2]), n) for n, r in rows.items())
    print(f"done: {len(rows)} values of n in {wall:.0f}s wall ({sum(times.values()):.0f}s CPU in workers)")
    print(f"  smallest min_M = {minM[0]:.6f} at n={minM[1]};  non-positive: "
          f"{[n for n, r in rows.items() if float(r[2]) <= 0]}")
    print(f"  switch points not strictly increasing in (1/2,1): {nonmono}")
    print(f"  max |balance residual| at the chosen p_m: {maxres:.2e}")
    print(f"  direct-summation cross-checks: {len(cross)}, max rel diff {max(r[4] for r in cross):.2e}")
