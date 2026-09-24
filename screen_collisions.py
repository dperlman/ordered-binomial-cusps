"""screen_collisions.py -- the accelerated collision search, licensed by Fact C.

TIER 2.  Exhaustive IF AND ONLY IF Fact C holds (RESEARCH_LOG, 2026-09-23).  Fact C says a tie
point whose radical SIMPLIFIES must satisfy  min(width, band) <= g(j) := j - (largest prime <= j),
because a prime in (max(i, n-i, j/2), j] divides C(n,i)/C(n,j) to exponent exactly 1, which forbids
simplification.  So instead of all ~n^2/4 pairs we examine only ~2*sum_j g(j) of them -- mean g(j)
is ~10 at j = 10^6 where n/4 would be 250,000.

Then Fact B (proved): any collision has at least one simplifying partner.  So for each simplifying
point we look for its partner, and that is the whole search.  The partner must share the TRUE ROOT
ORDER M (least M with rho^M rational) and the same reduced value.  If the partner does not itself
simplify its width IS M, and Fact A (proved: rho is strictly increasing in i at fixed width) makes
that a binary search.  Simplifying points are also compared against each other directly.

    .venv/bin/python screen_collisions.py --nmax 100000 --workers 8
    .venv/bin/python screen_collisions.py --nmax 400 --verify      # cross-check vs the full scan
"""
import argparse, sys, time
from math import gcd, lgamma
from multiprocessing import Pool, cpu_count
import numpy as np
from check_collisions import Valuation, primes_upto, digit_sum, digit_sum_vec

def prev_prime_array(n):
    isp = np.zeros(n+1, bool); isp[primes_upto(n)] = True
    return np.maximum.accumulate(np.where(isp, np.arange(n+1), 0))

N_STAGE1 = 8             # primes in the vectorised whole-array stage.  8 measured best:
                         # more primes cut survivors (2400 -> 246 at 24 primes) but the whole-array
                         # cost grows faster than the saving on the survivor stages.

def _vp_at(n, p, ks):
    """v_p(C(n,k)) for an ARRAY of k, by Legendre, without building the whole length-n array."""
    return ((digit_sum_vec(ks, p) + digit_sum_vec(n - ks, p) - digit_sum(n, p)) // (p - 1))

def _ranges_to_pairs(jv, lo, hi):
    """Vectorised: for each j, the integers i in [lo, hi] -> flat (i, j) arrays."""
    cnt = np.maximum(hi - lo + 1, 0)
    if cnt.sum() == 0:
        return np.empty(0, np.int64), np.empty(0, np.int64)
    jj = np.repeat(jv, cnt)
    start = np.cumsum(cnt) - cnt
    off = np.arange(cnt.sum()) - np.repeat(start, cnt)
    ii = np.repeat(lo, cnt) + off
    return ii, jj

def simplifying_at(n, val, prevp, restrict=True):
    """Every simplifying tie point of n.  restrict=True applies Fact C's g(j) window.

    Fully vectorised: the candidate list is built with np.repeat, and the d>1 test is done by
    gcd-ing the exponents of the first few primes.  d = gcd(m, gcd_p e_p) must divide
    gcd(m, gcd of ANY subset of the e_p), so a cheap gcd over 8 primes rejects almost everything
    before any per-candidate Python work.
    """
    jv = np.arange(2, n+1, dtype=np.int64)
    g = (jv - prevp[jv]) if restrict else jv.copy()
    one = np.ones_like(jv)
    loA = np.maximum(np.maximum(jv - g, n - jv + 1), one)      # width <= g
    hiA = jv - 1
    loB = np.maximum(n - jv + 1, one)                          # band  <= g
    hiB = np.minimum(n - jv + g, jv - 1)
    iA, jA = _ranges_to_pairs(jv, loA, hiA)
    iB, jB = _ranges_to_pairs(jv, loB, hiB)
    i = np.concatenate([iA, iB]); j = np.concatenate([jA, jB])
    if len(i) == 0: return []
    # NOTE: windows A and B can overlap, but the overlap is ~0.1% of candidates and testing one
    # twice is harmless, so we do NOT deduplicate here -- np.unique on the full array was 26% of
    # the per-n cost.  The few duplicates are removed from the (tiny) survivor list instead.
    m = j - i
    keep = (m >= 2) & (i + j > n) & (i >= 1) & (j <= n)
    i, j, m = i[keep], j[keep], m[keep]
    if len(i) == 0: return []
    # ESCALATING PREFILTER.  d = gcd(m, gcd over ALL primes of e_p).  The gcd over any SUBSET is a
    # multiple of the true gcd, so gcd(m, G_subset) == 1 is a VALID rejection.  Eight primes leaves
    # ~2^-8 of candidates alive, which is still ~80 per n and each of those costs a full Python
    # sweep over pi(n) primes -- that was 94% of the runtime.  Escalate in stages instead: every
    # stage is vectorised, and almost nothing reaches Python.
    idx = np.arange(len(i))
    G = np.zeros(len(i), np.int64)
    for lo_p, hi_p, whole in ((0, N_STAGE1, True), (N_STAGE1, 64, False), (64, 256, False)):
        ps = val.primes[lo_p:hi_p]
        if len(ps) == 0 or len(idx) == 0: break
        for p in ps:
            pi = int(p)
            if whole:                       # early stage: whole-array vec is cached and reused
                V = val.vec(pi)
                e = V[i[idx]] - V[j[idx]]
            else:                           # later stages: few survivors, so evaluate POINTWISE
                e = _vp_at(n, pi, i[idx]) - _vp_at(n, pi, j[idx])
            G[idx] = np.gcd(G[idx], e)
        idx = idx[np.gcd(m[idx], G[idx]) > 1]
    alive = np.zeros(len(i), bool); alive[idx] = True
    out = []
    seen = set()
    for t in np.flatnonzero(alive):                            # survivors only: full prime sweep
        ii, jj_, mm = int(i[t]), int(j[t]), int(m[t])
        e = [val.v(int(p), ii) - val.v(int(p), jj_) for p in val.primes]
        ge = 0
        for x in e: ge = gcd(ge, x)
        d = gcd(mm, ge)
        if d > 1 and (ii, jj_) not in seen:
            seen.add((ii, jj_))
            out.append((ii, jj_, mm, d, mm//d, tuple(x//d for x in e)))
    return out

def find_partner(n, val, lnC, i, j, m, M, target):
    """Is there another tie point with the same true root order M and the same reduced value?"""
    hits = []
    lnrho = (lnC[i] - lnC[j]) / m
    lo, hi = 1, n - M                                   # width-M pairs (i2, i2+M), i2+i2+M > n
    lo = max(lo, (n - M)//2)
    if hi >= lo:
        a, b = lo, hi                                   # binary search: rho increases with i2
        while a < b:
            mid = (a + b)//2
            if (lnC[mid] - lnC[mid+M])/M < lnrho: a = mid + 1
            else: b = mid
        for i2 in range(max(lo, a-4), min(hi, a+4)+1):  # exact-verify a small window
            j2 = i2 + M
            if (i2, j2) == (i, j) or i2 + j2 <= n: continue
            if all(val.v(int(p), i2) - val.v(int(p), j2) == t
                   for p, t in zip(val.primes, target)):
                hits.append((i2, j2))
    return hits

def one_n(args):
    n, verify = args[0], args[1]
    val = Valuation(n)
    prevp = prev_prime_array(n)
    lnC = np.array([lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1) for k in range(n+1)])
    S = simplifying_at(n, val, prevp)
    if verify:
        full = simplifying_at(n, val, prevp, restrict=False)
        if {(x[0], x[1]) for x in S} != {(x[0], x[1]) for x in full}:
            return n, len(S), [("FACT-C-MISS", sorted({(x[0],x[1]) for x in full} -
                                                      {(x[0],x[1]) for x in S}))]
    hits = []
    for (i, j, m, d, M, target) in S:                   # partner of each simplifying point
        for (i2, j2) in find_partner(n, val, lnC, i, j, m, M, target):
            hits.append((n, i, j, i2, j2))
    for a in range(len(S)):                             # simplifying vs simplifying
        for b in range(a+1, len(S)):
            if S[a][4] == S[b][4] and S[a][5] == S[b][5]:
                hits.append((n, S[a][0], S[a][1], S[b][0], S[b][1]))
    return n, len(S), hits, [(n, x[0], x[1], x[2], x[0]+x[1]-n, x[4]) for x in S]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=1000)
    ap.add_argument("--nmin", type=int, default=3)
    ap.add_argument("--workers", default="auto")
    ap.add_argument("--save", metavar="CSV",
                    help="write the catalogue of simplifying tie points to this CSV")
    ap.add_argument("--verify", action="store_true",
                    help="also run the UNRESTRICTED scan at each n and confirm Fact C missed nothing")
    a = ap.parse_args()
    w = cpu_count() if a.workers == "auto" else int(a.workers)
    ns = list(range(a.nmax, a.nmin-1, -1))
    print(f"screening n={a.nmin}..{a.nmax} on {w} workers"
          + ("  [--verify: Fact C checked against the full scan at every n]" if a.verify else ""),
          flush=True)
    t0 = time.time(); tot_s = 0; allhits = []; done = 0
    save = [] if a.save else None
    with Pool(w) as pool:
        for n, ns_count, hits, rows in pool.imap_unordered(one_n, [(n, a.verify) for n in ns],
                                                     chunksize=4):
            tot_s += ns_count; done += 1; allhits += hits
            if save is not None: save.extend(rows)
            if hits or done % 2000 == 0:
                print(f"  [{done}/{len(ns)}  {time.time()-t0:.0f}s]  n={n}  "
                      f"simplifying so far {tot_s}" + (f"  HITS {hits}" if hits else ""), flush=True)
    print(f"\n{tot_s} simplifying tie points found over {len(ns)} values of n "
          f"in {time.time()-t0:.0f}s")
    miss = [h for h in allhits if h and h[0] == "FACT-C-MISS"]
    real = [h for h in allhits if not (h and h[0] == "FACT-C-MISS")]
    if miss:
        print(f"\n*** FACT C FAILED at {len(miss)} values of n ***")
        for m_ in miss[:10]: print("   ", m_)
    if save is not None:
        save.sort()
        with open(a.save, "w") as fh:
            fh.write("n,i,j,width,band,true_root_order\n")
            for r in save: fh.write(",".join(map(str, r)) + "\n")
        print(f"catalogue of {len(save)} simplifying tie points -> {a.save}")
    print("\n" + "="*70)
    if real:
        print(f"*** COLLISIONS FOUND: {len(real)} ***")
        for h in real[:50]: print("   ", h)
        sys.exit(1)
    print("NO COLLISIONS found by the Fact C screen.")
    if miss: sys.exit(1)

if __name__ == "__main__":
    main()
