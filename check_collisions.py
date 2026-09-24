"""check_collisions.py -- are there two DIFFERENT (i,j) pairs with the same tie point p*?

A tie point is p* = rho/(1+rho) with rho^(j-i) = C(n,i)/C(n,j), so two pairs (i,j) and (k,l)
collide exactly when

    (C(n,i)/C(n,j))^(l-k)  ==  (C(n,k)/C(n,l))^(j-i)                                       (*)

That is an identity between integers, so it can be settled EXACTLY -- and cheaply, without ever
forming the numbers.  Take the p-adic valuation of both sides: (*) holds iff, for every prime
p <= n,

    (l-k) * (v_p C(n,i) - v_p C(n,j))  ==  (j-i) * (v_p C(n,k) - v_p C(n,l))

with v_p C(n,k) = (S_p(k) + S_p(n-k) - S_p(n)) / (p-1)   [Legendre],  S_p = digit sum in base p.
Almost every candidate fails at one of the first few primes, so the early exit makes this fast.

Double precision CANNOT settle this on its own at large n: the closest distinct tie points at
n=6000 are 4.4e-16 apart, which is 4 ulps, and the computed p* itself carries an error of order
q*3*eps*max(lnC)/(j-i) -- up to ~4e-12 for narrow pairs.  So the numeric pass is only a SCREEN: it
finds every pair closer than --tol, and each of those is then decided exactly.

    .venv/bin/python check_collisions.py --parquet          # the complete tie-point dumps
    .venv/bin/python check_collisions.py --cusps            # cusp-cusp only, every n <= 5000
    .venv/bin/python check_collisions.py --exhaustive 2000  # EVERY tie point of every n <= 2000
"""
import argparse, glob, os, sys
from multiprocessing import Pool, cpu_count
import numpy as np

def primes_upto(n):
    s = np.ones(n+1, bool); s[:2] = False
    for p in range(2, int(n**0.5)+1):
        if s[p]: s[p*p::p] = False
    return np.flatnonzero(s)

def digit_sum(x, p):
    s = 0
    while x: s += x % p; x //= p
    return s

def digit_sum_vec(x, p):
    """Digit sum in base p, elementwise, for an integer array."""
    x = x.astype(np.int64, copy=True); s = np.zeros_like(x)
    while x.any():
        s += x % p; x //= p
    return s

class Valuation:
    """v_p(C(n,k)) for one n, any prime p and any k, via Legendre's formula.

    vec(p) gives the whole array k=0..n at once, which is what makes the candidate prefilter
    vectorised: most candidates are separated by the very first prime, so testing 2,3,5,7,...
    in numpy and only sending survivors to the scalar sweep turns an O(candidates * primes)
    Python loop into a handful of array ops.
    """
    def __init__(self, n):
        self.n = n; self.primes = primes_upto(n)
        self._Sn = {}; self._vec = {}
    def v(self, p, k):
        n = self.n
        Sn = self._Sn.get(p)
        if Sn is None: Sn = self._Sn[p] = digit_sum(n, p)
        return (digit_sum(k, p) + digit_sum(n-k, p) - Sn) // (p-1)
    def vec(self, p):
        a = self._vec.get(p)
        if a is None:
            k = np.arange(self.n+1)
            a = self._vec[p] = ((digit_sum_vec(k, p) + digit_sum_vec(self.n-k, p)
                                 - digit_sum(self.n, p)) // (p-1))
        return a

def same_tie_point(val, i, j, k, l):
    """EXACT: do (i,j) and (k,l) have the same p*?  (*) above, prime by prime."""
    a, b = l - k, j - i
    for p in val.primes:
        if a*(val.v(p, i) - val.v(p, j)) != b*(val.v(p, k) - val.v(p, l)):
            return False
    return True

def screen(n, i, j, p, tol, val=None, label=""):
    """Sort by p*, exact-check EVERY pair closer than tol.  Returns (collisions, min gap, #pairs).

    Every pair, not just sorted neighbours.  Until 2026-09-24 this compared only consecutive
    pairs, which is NOT sound: two tie points A and C that collide exactly still have computed p*
    differing by rounding, and a near-miss B can land between them in float order; then (A,B) and
    (B,C) are tested and (A,C) never is.  (Caught by external review.)  The fix walks offsets
    k = 1, 2, ... and takes every pair (a, a+k) with p[a+k] - p[a] < tol.  It stops at the first
    k with no such pair, which is safe because p is sorted: p[a+k+1]-p[a] >= p[a+k]-p[a] >= tol.
    """
    o = np.argsort(p, kind='stable'); i, j, p = i[o], j[o], p[o]
    g = np.diff(p)
    ia, ib = [], []
    k = 1
    while k < len(p):
        idx = np.flatnonzero(p[k:] - p[:-k] < tol)
        if not len(idx): break
        ia.append(idx); ib.append(idx + k)
        k += 1
    hits = []
    if not ia:
        return hits, (g.min() if len(g) else np.inf), 0
    a_ = np.concatenate(ia); b_ = np.concatenate(ib)
    n_cand = len(a_)
    val = val or Valuation(n)
    i1, j1 = i[a_].astype(np.int64), j[a_].astype(np.int64)
    i2, j2 = i[b_].astype(np.int64), j[b_].astype(np.int64)
    A, B = (j2 - i2), (j1 - i1)                     # (l-k) and (j-i)
    alive = np.ones(n_cand, bool)
    for q in val.primes[:12]:                       # vectorised prefilter on the small primes
        V = val.vec(int(q))
        alive &= (A*(V[i1] - V[j1]) == B*(V[i2] - V[j2]))
        if not alive.any(): break
    for c in np.flatnonzero(alive):                 # survivors: full prime sweep, exact
        if same_tie_point(val, int(i1[c]), int(j1[c]), int(i2[c]), int(j2[c])):
            hits.append((n, int(i1[c]), int(j1[c]), int(i2[c]), int(j2[c]), float(p[a_[c]])))
    return hits, (g.min() if len(g) else np.inf), n_cand

def tie_points(n, lnC):
    """Every tie point of n: i<j<=n, i+j>n.  Returns (i, j, p*)."""
    lo = np.maximum(np.arange(1, n) + 1, n - np.arange(1, n) + 1)
    cnt = np.maximum(0, n - lo + 1)
    ii = np.repeat(np.arange(1, n), cnt)
    jj = np.concatenate([np.arange(lo[t], n+1) for t in range(n-1) if cnt[t] > 0])
    lnr = (lnC[ii] - lnC[jj]) / (jj - ii)
    return ii, jj, 1.0/(1.0 + np.exp(-lnr))

def _one_n(args):
    """One n, for the worker pool: generate every tie point, screen, exact-check."""
    from math import lgamma
    n, tol = args
    lnC = np.array([lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1) for k in range(n+1)])
    i, j, p = tie_points(n, lnC)
    hits, gmin, ncand = screen(n, i, j, p, tol)
    return n, len(i), gmin, ncand, hits

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parquet", action="store_true")
    ap.add_argument("--cusps", action="store_true")
    ap.add_argument("--exhaustive", type=int, metavar="NMAX")
    ap.add_argument("--data", default="data/ties")
    ap.add_argument("--cusps-csv", default="cusps/cusps_all.csv")
    ap.add_argument("--workers", default="auto")
    ap.add_argument("--tol", type=float, default=1e-9,
                    help="exact-check every consecutive pair closer than this (default 1e-9, "
                         "~250x the worst numeric error in p*)")
    a = ap.parse_args()
    total_hits = []

    if a.parquet:
        import pyarrow.parquet as pq
        print(f"{'n':>6} {'tie points':>12} {'min gap':>12} {'candidates':>11} {'collisions':>11}")
        for d in sorted(glob.glob(os.path.join(a.data, 'n=*'))):
            f = os.path.join(d, 'part.parquet')
            if not os.path.exists(f): continue
            n = int(os.path.basename(d).split('=')[1])
            t = pq.read_table(f, columns=['i', 'j', 'pstar'])
            i = t['i'].to_numpy(); j = t['j'].to_numpy(); p = t['pstar'].to_numpy()
            keep = i < j                       # drop the p=1/2 axis row (i=0, j=n)
            keep &= (i + j) > n
            hits, gmin, ncand = screen(n, i[keep], j[keep], p[keep], a.tol)
            total_hits += hits
            print(f"{n:6d} {int(keep.sum()):12,d} {gmin:12.3e} {ncand:11d} {len(hits):11d}",
                  flush=True)

    if a.cusps:
        import csv
        print("\nCUSP-CUSP collisions across every n in", a.cusps_csv)
        cur = None; I = []; J = []; P = []; nn = 0; mg = np.inf; nc = 0
        def flush(nv):
            nonlocal mg, nc
            if nv is None or len(I) < 2: return []
            h, g, c = screen(nv, np.array(I), np.array(J), np.array(P), a.tol)
            mg = min(mg, g); nc += c
            return h
        with open(a.cusps_csv) as fh:
            rd = csv.reader(fh); next(rd)
            for r in rd:
                n = int(r[0])
                if n != cur:
                    total_hits += flush(cur); nn += 1
                    cur = n; I, J, P = [], [], []
                I.append(int(r[1])); J.append(int(r[2])); P.append(float(r[3]))
        total_hits += flush(cur)
        print(f"  {nn} values of n scanned; min cusp-cusp gap {mg:.3e}; "
              f"{nc} candidates; {len(total_hits)} collisions")

    if a.exhaustive:
        import time
        w = cpu_count() if a.workers == "auto" else int(a.workers)
        print(f"\nEVERY tie point of every n <= {a.exhaustive}   ({w} workers)")
        print(f"{'n':>6} {'tie points':>12} {'min gap':>12} {'candidates':>11} {'collisions':>11}")
        ns = list(range(a.exhaustive, 2, -1))          # largest n first: better load balance
        grand = 0; done = 0; t0 = time.time(); worst = {}
        with Pool(w) as pool:
            for n, cnt, gmin, ncand, hits in pool.imap_unordered(_one_n,
                    [(n, a.tol) for n in ns], chunksize=4):
                grand += cnt; done += 1; total_hits += hits
                worst[n] = (cnt, gmin, ncand)
                if hits or done % 500 == 0:
                    el = time.time() - t0
                    print(f"{n:6d} {cnt:12,d} {gmin:12.3e} {ncand:11d} {len(hits):11d}"
                          f"   [{done}/{len(ns)}, {el:.0f}s]", flush=True)
        for n in sorted(worst)[::max(1, len(worst)//12)]:
            cnt, gmin, ncand = worst[n]
            print(f"{n:6d} {cnt:12,d} {gmin:12.3e} {ncand:11d}")
        print(f"  {grand:,} tie points checked in total in {time.time()-t0:.0f}s")

    print("\n" + "="*70)
    if total_hits:
        print(f"COLLISIONS FOUND: {len(total_hits)}")
        for h in total_hits[:50]: print("   ", h)
        sys.exit(1)
    print("NO COLLISIONS: every pair that double precision could not separate was checked")
    print("exactly (prime-by-prime) and proved distinct.")

if __name__ == "__main__":
    main()
