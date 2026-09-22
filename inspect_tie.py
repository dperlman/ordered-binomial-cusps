"""inspect_tie.py -- dissect one tie point's ranking and its re-ranking clusters.

Shows, for a given (n,i,j), what the sharpened CHECK trigger sees: the masses in RANK order
(increasing), the relative gap between each pair of neighbours, the resolution threshold from
binom_core._err_bounds, and which adjacencies double precision cannot resolve.  Maximal runs of
unresolvable neighbours are the "clusters"; within a cluster of size c any rank can move by c-1,
so the cluster contributes (c-1)*sum_{k in C} f_k|k - n p*| to the bound that is added to MARGIN.

Run: .venv/bin/python inspect_tie.py 2590 791 2243
     .venv/bin/python inspect_tie.py 2590 791 2243 --top 40 --csv ranking.csv
"""
import argparse
from math import log
import numpy as np
import binom_core as core

def dissect(n, i, j):
    lnC = core.lnC_arr(n)
    p = 1.0/(1.0 + np.exp(-(lnC[i]-lnC[j])/(j-i))); q = 1.0 - p
    lnp, lnq = log(p), log(q)
    md = min(int((n+1)*p), n)
    k = np.arange(n+1)
    f = np.exp(lnC + k*lnp + (n-k)*lnq)
    f[f < core.TINY] = 0.0
    f[j] = f[i]                                   # the tie is exact
    f = f/f.sum()                                 # same normalisation as the kernel
    dp, df0, dstep = core._err_bounds(n, lnC, i, j, md, q, lnp, lnq)
    order = np.lexsort((k, f))                    # increasing mass, index as tiebreak
    a, b = order[:-1], order[1:]
    gap = np.where(f[b] > 0, (f[b]-f[a])/np.where(f[b] > 0, f[b], 1.0), np.inf)
    thr = 2*df0 + dstep*(np.abs(b-md) + np.abs(a-md)) + dp*np.abs(b-a)/q
    return dict(p=p, q=q, md=md, f=f, order=order, gap=gap, thr=thr,
                dp=dp, df0=df0, dstep=dstep, lnC=lnC)

def clusters(d, i, j):
    order, gap, thr, f, p = d['order'], d['gap'], d['thr'], d['f'], d['p']
    out = []; cur = [int(order[0])]
    for t in range(len(order)-1):
        if gap[t] <= thr[t]: cur.append(int(order[t+1]))
        else: out.append(cur); cur = [int(order[t+1])]
    out.append(cur)
    return [c for c in out if len(c) > 1]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('n', type=int); ap.add_argument('i', type=int); ap.add_argument('j', type=int)
    ap.add_argument('--top', type=int, default=25, help='how many of the largest masses to list')
    ap.add_argument('--csv', help='also write the whole ranking to this CSV')
    a = ap.parse_args()
    n, i, j = a.n, a.i, a.j
    d = dissect(n, i, j); f, order, p, q, md = d['f'], d['order'], d['p'], d['q'], d['md']
    N = n + 1
    print(f"n={n}  i={i}  j={j}   p*={p:.15f}   n*p*={n*p:.2f}   mode={md}   "
          f"sqrt(npq)={np.sqrt(n*p*q):.2f}")
    print(f"tied pair f(i)=f(j) = {f[i]:.6e}   kappa=(j-i)f(i) = {(j-i)*f[i]:.6e}")
    print(f"largest mass f({int(f.argmax())}) = {f.max():.6e}")
    print(f"error bounds: dp={d['dp']:.3e}  df0={d['df0']:.3e}  dstep={d['dstep']:.3e}")
    print(f"\nRanking, largest {a.top} masses (rank = position in INCREASING order, {N} masses).")
    print("'gap' is the relative gap down to the next mass below; unresolvable when gap <= thr.\n")
    print(f"{'rank':>6} {'k':>6} {'f(k)':>25} {'gap':>12} {'thr':>11}  {'':1}")
    for pos in range(N-1, max(N-1-a.top, 0), -1):
        t = pos - 1
        mark = '  <== UNRESOLVABLE' if d['gap'][t] <= d['thr'][t] else ''
        print(f"{pos:>6} {int(order[pos]):>6} {f[order[pos]]:>25.17e} "
              f"{d['gap'][t]:>12.3e} {d['thr'][t]:>11.3e}{mark}")
    cl = clusters(d, i, j)
    print(f"\n{len(cl)} cluster(s) of size > 1:")
    tot = 0.0
    for c in cl:
        if set(c) <= {i, j}:
            print(f"  {{i,j}} = {{{i},{j}}}: exact tie, left-limit ranking is by DEFINITION "
                  f"-> contributes 0"); continue
        if max(f[x] for x in c) == 0.0:
            print(f"  zero block: {len(c)} masses below TINY -> contributes 0"); continue
        s = sum(f[x]*abs(x - n*p) for x in c); b = (len(c)-1)*s; tot += b
        print(f"  size {len(c)}: k = {c}")
        for x in c:
            print(f"      k={x:6d}  f={f[x]:.17e}  k-n*p*={x-n*p:+9.2f}  "
                  f"{100*f[x]/f.max():5.1f}% of peak")
        print(f"      contributes (c-1)*sum f|a| = {b:.6e}")
    *_, Sm, kappa, F3, tag_old, tag_new, rbnd = core._one_tie(
        n, d['lnC'], i, j, np.empty(n+1), np.empty(n+1, np.int64))
    nm = {core.TAG_NOT: 'NOT', core.TAG_MIN: 'MIN', core.TAG_CHECK: 'CHECK'}
    print(f"\nre-ranking bound  = {rbnd:.6e}   (kernel: {rbnd:.6e})")
    print(f"S_- = {Sm:.9e}   S_+ = {Sm+kappa:.9e}   MARGIN = {core.MARGIN:.0e}")
    print(f"tag: legacy GAP rule -> {nm[tag_old]},  sharpened -> {nm[tag_new]}")
    if a.csv:
        with open(a.csv, 'w') as fh:
            fh.write("rank,k,f,gap_to_next_below,threshold,unresolvable\n")
            for pos in range(N-1, 0, -1):
                t = pos - 1
                fh.write(f"{pos},{int(order[pos])},{f[order[pos]]:.17e},{d['gap'][t]:.6e},"
                         f"{d['thr'][t]:.6e},{int(d['gap'][t] <= d['thr'][t])}\n")
        print(f"\nfull ranking written to {a.csv}")

if __name__ == '__main__':
    main()
