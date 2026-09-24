"""verify_fact_c.py -- independent check that the Fact C window loses no reducing tie point.

For each n it finds EVERY reducing ("simplifying") tie point by brute force -- every width m from
2 to n-1, every valid i -- with no use of Fact C at all, and compares the result against
  (a) screen_collisions.simplifying_at(restrict=True), i.e. what the Fact C screen examines, and
  (b) the rows for that n in the saved catalogue (--catalogue), if given.
Any reducing point found here but missing from (a) or (b) is a failure of Fact C or of its
implementation.

Deliberately independent of screen_collisions.py's window logic.  It shares only Legendre's
formula (check_collisions.Valuation), which both rest on.  Organised by WIDTH so memory is O(n)
per n rather than the O(n^2) of screen_collisions --verify, so it can run at n ~ 10^5.

Test for "reducing": d = gcd(m, gcd_p e_p) > 1, where e_p = v_p C(n,i) - v_p C(n,j).  Equivalently
some prime q | m divides every e_p.  For each prime q | m the candidates are prefiltered on the
first few primes (compressing after each), then survivors get the full sweep over all p <= n.

    .venv/bin/python verify_fact_c.py --nmin 39950 --nmax 40050 --workers 8 \
        --catalogue simplifying_n100000.csv
"""
import argparse, sys, time
from math import gcd
from multiprocessing import Pool, cpu_count
import numpy as np
from check_collisions import Valuation

def spf_sieve(n):
    s = np.arange(n+1)
    for p in range(2, int(n**0.5)+1):
        if s[p] == p:
            blk = s[p*p::p]
            blk[blk == np.arange(p*p, n+1, p)] = p
    return s

def prime_divisors(m, spf):
    out = []
    while m > 1:
        q = int(spf[m]); out.append(q)
        while m % q == 0: m //= q
    return out

def brute_force(n, val, npre=6):
    """Every reducing tie point of n, with no appeal to Fact C."""
    V = [val.vec(int(p)) for p in val.primes[:npre]]
    spf = spf_sieve(n)
    found = set()
    for m in range(2, n):
        lo, hi = (n - m)//2 + 1, n - m
        if hi < lo: continue
        i = np.arange(lo, hi+1); j = i + m
        for q in prime_divisors(m, spf):
            ia, ja = i, j
            for Vp in V:
                keep = (Vp[ia] - Vp[ja]) % q == 0
                ia, ja = ia[keep], ja[keep]
                if not len(ia): break
            for t in range(len(ia)):
                ii, jj = int(ia[t]), int(ja[t])
                if all((val.v(int(p), ii) - val.v(int(p), jj)) % q == 0 for p in val.primes):
                    found.add((ii, jj))
    return found

def one_n(n):
    import screen_collisions as sc
    t0 = time.time()
    val = Valuation(n)
    full = brute_force(n, val)
    restricted = {(x[0], x[1]) for x in sc.simplifying_at(n, val, sc.prev_prime_array(n))}
    return n, full, restricted, time.time() - t0

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int)
    ap.add_argument("--nmax", type=int)
    ap.add_argument("--ns", type=int, nargs="*", default=[],
                    help="extra individual n to check, e.g. where n+1 is highly composite")
    ap.add_argument("--workers", default="auto")
    ap.add_argument("--catalogue", help="CSV from screen_collisions --save, to compare against")
    a = ap.parse_args()
    want = set(a.ns)
    if a.nmin is not None and a.nmax is not None: want |= set(range(a.nmin, a.nmax + 1))
    cat = {}
    if a.catalogue:
        for line in open(a.catalogue).read().splitlines()[1:]:
            n, i, j = map(int, line.split(',')[:3])
            if n in want: cat.setdefault(n, set()).add((i, j))
    w = cpu_count() if a.workers == "auto" else int(a.workers)
    ns = sorted(want, reverse=True)                     # largest first: better load balance
    print(f"brute-force verification of the Fact C window at {len(ns)} values of n "
          f"({min(ns)}..{max(ns)}), {w} workers", flush=True)
    t0 = time.time(); bad = []; tot = 0
    with Pool(w) as pool:
        for n, full, restricted, dt in pool.imap_unordered(one_n, ns):
            tot += len(full)
            miss_screen = full - restricted
            extra_screen = restricted - full
            line = f"  n={n}: {len(full)} reducing (brute force), {dt:.0f}s"
            if a.catalogue:
                c = cat.get(n, set())
                line += f", catalogue {len(c)}"
                if c != full: bad.append((n, "catalogue", sorted(full - c), sorted(c - full)))
            if miss_screen or extra_screen:
                bad.append((n, "screen", sorted(miss_screen), sorted(extra_screen)))
                line += "  *** MISMATCH ***"
            print(line, flush=True)
    print(f"\n{len(ns)} values of n, {tot} reducing tie points by brute force, {time.time()-t0:.0f}s")
    print("=" * 70)
    if bad:
        print("*** MISMATCHES -- Fact C or its implementation is wrong here ***")
        for b in bad: print("   ", b)
        sys.exit(1)
    print("AGREE at every n: brute force, the Fact C screen" +
          (" and the catalogue" if a.catalogue else "") + " find exactly the same reducing points.")

if __name__ == "__main__":
    main()
