"""
analyze_cusps.py -- analyses over a merged cusp table (cusps_all.csv or slim parts).

    python analyze_cusps.py cusps/cusps_all.csv [--nmax N] [--out analysis/]

Writes to --out:
  per_n_summary.csv      one row per n: counts, extremes, gap to E(1/2)
  cusps_F3_negative.csv  every F3<0 cusp with nearest-cusp distance metrics
  all_neighbors.csv      nearest-cusp distance metrics for ALL cusps (for baselines)
  summary.txt            headline numbers and outlier lists
Needs only numpy.  Tie points are regenerated per n (vectorised), not read from disk.
"""
import argparse, csv, os, sys, time
import numpy as np
from math import lgamma
from collections import defaultdict

from binom_core import lnC_arr, E_half as _E_half     # single implementation, see binom_core.py

def tie_points(n, lnC):
    """all (p*, i, j) with 0<=i<j<=n, i+j>n, sorted by p*

    j<=n matters here: the pairs (i,n) are real tie points and they sit between the others in p*
    order, so excluding them would corrupt every nearest-neighbour gap and intervening-tie count.
    """
    i, j = np.triu_indices(n+1, 1)
    m = (i >= 1) & (i+j > n); i, j = i[m], j[m]
    p = 1/(1+np.exp(-(lnC[i]-lnC[j])/(j-i)))
    o = np.argsort(p, kind='stable'); return p[o], i[o], j[o]

def E_half(n, lnC=None):
    """E(n,1/2) from binom_core (normalised masses).  lnC is accepted and ignored, for callers
    that still pass it."""
    return _E_half(n)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", nargs='+'); ap.add_argument("--nmax", type=int, default=10**9)
    ap.add_argument("--out", default="analysis"); a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rows = []
    for fn in a.csv:
        with open(fn) as f:
            for r in csv.DictReader(f):
                if r['pstar'] == '' or int(r['n']) > a.nmax: continue
                rows.append(r)
    byn = defaultdict(list)
    for r in rows: byn[int(r['n'])].append(r)
    print(f"{len(rows)} cusps over {len(byn)} values of n (max n={max(byn)})", flush=True)
    t0 = time.time()
    per_n = []; nb_rows = []
    for n in sorted(byn):
        lnC = lnC_arr(n); p, ti, tj = tie_points(n, lnC)
        idx = {(int(x), int(y)): t for t, (x, y) in enumerate(zip(ti, tj))}
        C = byn[n]; cidx = sorted((idx[(int(r['i']), int(r['j']))], r) for r in C)
        Eh = E_half(n, lnC)
        Es = np.array([float(r['E']) for r in C]); F3s = np.array([float(r['F3']) for r in C])
        ps = np.array([float(r['pstar']) for r in C])
        jumps = np.array([float(r['slope_right'])-float(r['slope_left']) for r in C])
        # nearest-cusp metrics
        gaps = []
        info = []
        for a_, (t, r) in enumerate(cidx):
            cand = []
            if a_ > 0:            tp, rp = cidx[a_-1]; cand.append((p[t]-p[tp], t-tp-1, rp))
            if a_+1 < len(cidx):  tn, rn = cidx[a_+1]; cand.append((p[tn]-p[t], tn-t-1, rn))
            if not cand: continue
            g, between, nb = min(cand, key=lambda c: c[0]); gaps.append(g); info.append((r, g, between, nb))
        gaps = np.array(gaps); med = np.median(gaps) if len(gaps) else float('nan')
        for (r, g, between, nb) in info:
            pct = 100*np.mean(gaps <= g)
            nb_rows.append(dict(n=n, i=int(r['i']), j=int(r['j']), pstar=float(r['pstar']), F3=float(r['F3']),
                                F3_sign=r['F3_sign'], nb_i=int(nb['i']), nb_j=int(nb['j']), nb_width=int(nb['j'])-int(nb['i']),
                                nb_F3_sign=nb['F3_sign'], gap=g, n_gap=n*g, intervening_ties=between,
                                percentile_at_n=pct, gap_over_median=g/med if med > 0 else float('nan'),
                                certified_by=r.get('certified_by', '')))
        neg = F3s < 0
        per_n.append(dict(n=n, n_cusps=len(C), n_negF3=int(neg.sum()), max_pstar=ps.max(), min_pstar=ps.min(),
                          min_F3=F3s.min(), E_half=Eh, min_E_minus_Ehalf=(Es-Eh).min(),
                          argmin_i=int(C[int(np.argmin(Es))]['i']), argmin_j=int(C[int(np.argmin(Es))]['j']),
                          median_slope_jump=float(np.median(jumps)), min_slope_jump=float(jumps.min()),
                          n_double_cusps=int(np.sum([b == 0 for (_, _, b, _) in info])),
                          n_iv=sum(1 for r in C if r.get('certified_by', '').startswith('iv'))))
        if n % 200 == 0: print(f"  n={n} done ({time.time()-t0:.0f}s)", flush=True)
    # write outputs
    with open(os.path.join(a.out, "per_n_summary.csv"), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(per_n[0].keys()), lineterminator="\n")
        w.writeheader(); w.writerows(per_n)
    for name, sel in (("all_neighbors.csv", nb_rows), ("cusps_F3_negative.csv", [r for r in nb_rows if r['F3_sign'] == '-'])):
        with open(os.path.join(a.out, name), 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(nb_rows[0].keys()), lineterminator="\n")
            w.writeheader(); w.writerows(sel)
    # summary
    neg = [r for r in nb_rows if r['F3_sign'] == '-']; pos = [r for r in nb_rows if r['F3_sign'] == '+']
    L = []
    L.append(f"cusps: {len(rows)}; F3<0: {len(neg)} ({100*len(neg)/len(rows):.2f}%); n range {min(byn)}..{max(byn)}")
    L.append(f"max cusp p* overall: {max(r['max_pstar'] for r in per_n):.5f}; by n (every 250): " +
             ", ".join(f"{r['n']}:{r['max_pstar']:.4f}" for r in per_n if r['n'] % 250 == 0))
    L.append(f"first n with F3<0 cusp: {next((r['n'] for r in per_n if r['n_negF3'] > 0), None)}")
    L.append(f"F3<0 fraction by n-range: " + ", ".join(
        f"{lo}-{hi}: {sum(r['n_negF3'] for r in per_n if lo<=r['n']<=hi)}/{sum(r['n_cusps'] for r in per_n if lo<=r['n']<=hi)}"
        for lo, hi in ((3,500),(501,1000),(1001,1500),(1501,2000),(2001,3000),(3001,5000)) if any(lo<=r['n']<=hi for r in per_n)))
    def q(v, x): return float(np.percentile(v, x)) if len(v) else float('nan')
    for name, S in (("F3<0", neg), ("F3>0", pos)):
        ng = np.array([r['n_gap'] for r in S]); bt = np.array([r['intervening_ties'] for r in S]); pc = np.array([r['percentile_at_n'] for r in S])
        rt = np.array([r['gap_over_median'] for r in S])
        L.append(f"{name} cusps (N={len(S)}): n*gap median {np.median(ng):.4f} [q90 {q(ng,90):.4f}, max {ng.max():.4f}]; "
                 f"intervening ties median {np.median(bt):.0f} [q90 {q(bt,90):.0f}, max {bt.max()}]; "
                 f"percentile-at-n median {np.median(pc):.0f} [q90 {q(pc,90):.0f}, max {pc.max():.0f}]; "
                 f"gap/median median {np.median(rt):.3f} [q90 {q(rt,90):.3f}, max {rt.max():.3f}]")
    L.append(f"F3<0 cusps that are double cusps (0 intervening ties): {sum(1 for r in neg if r['intervening_ties']==0)} / {len(neg)}"
             f"   (F3>0: {sum(1 for r in pos if r['intervening_ties']==0)} / {len(pos)})")
    L.append(f"nearest cusp of F3<0 cusps: width j-i median {np.median([r['nb_width'] for r in neg]):.0f}, max {max(r['nb_width'] for r in neg)}; "
             f"nearest has F3<0 itself: {sum(1 for r in neg if r['nb_F3_sign']=='-')}")
    neg_sorted = sorted(neg, key=lambda r: -r['percentile_at_n'])[:15]
    L.append("\nF3<0 cusps least close to another cusp (top 15 by percentile at same n):")
    L.append("n,i,j,pstar,F3,nearest(i,j),n*gap,intervening,percentile,gap/median")
    for r in neg_sorted:
        L.append(f"{r['n']},{r['i']},{r['j']},{r['pstar']:.6f},{r['F3']:.3f},({r['nb_i']},{r['nb_j']}),{r['n_gap']:.4f},{r['intervening_ties']},{r['percentile_at_n']:.0f},{r['gap_over_median']:.3f}")
    iv = [r for r in nb_rows if r['certified_by'].startswith('iv')]
    L.append(f"\ninterval-certified cusps: {len(iv)}; with F3<0: {sum(1 for r in iv if r['F3_sign']=='-')}")
    L.append(f"closest cusp to E(1/2): min over n of (E-E(1/2))*n = {min(r['min_E_minus_Ehalf']*r['n'] for r in per_n):.4f}; "
             f"always in first band (i+j=n+1)? {all(r['argmin_i']+r['argmin_j']==r['n']+1 for r in per_n)}")
    txt = "\n".join(L); print(txt)
    open(os.path.join(a.out, "summary.txt"), 'w').write(txt + "\n")

if __name__ == "__main__": main()
