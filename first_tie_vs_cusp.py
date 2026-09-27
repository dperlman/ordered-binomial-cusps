"""
first_tie_vs_cusp.py -- is the first tie point above p=1/2 always the first cusp above 1/2?

    .venv/bin/python first_tie_vs_cusp.py [--nmax 5000] [--workers 8] [--out analysis/]

"First" means the smallest p* > 1/2 (the degenerate tie at p=1/2 itself, where every mirror pair
ties, is excluded).  For each n the first tie point is found by an exhaustive scan of all pairs
0<=i<j<=n, i+j>n: p* is increasing in ln rho = (lnC(i)-lnC(j))/(j-i), so the scan minimises ln rho.
It also records the runner-up, and the gap between the two, so a float near-tie cannot silently pick
the wrong pair.  The first cusp is the smallest-p* row of cusps/cusps_all.csv (certified).

Writes analysis/first_tie_vs_cusp.csv, one row per n.
"""
import argparse, csv, os
import numpy as np
from numba import njit
from binom_core import lnC_arr

@njit(cache=True)
def _two_lowest(n, lnC):
    """(i, j, lnrho) of the smallest and second-smallest ln rho over tie points with i+j>n."""
    b1 = np.inf; b2 = np.inf; i1 = j1 = i2 = j2 = -1
    for i in range(1, n):
        for j in range(max(i + 1, n + 1 - i), n + 1):
            v = (lnC[i] - lnC[j])/(j - i)
            if v < b1:
                b2, i2, j2 = b1, i1, j1
                b1, i1, j1 = v, i, j
            elif v < b2:
                b2, i2, j2 = v, i, j
    return i1, j1, b1, i2, j2, b2

def first_tie(n):
    i1, j1, b1, i2, j2, b2 = _two_lowest(n, lnC_arr(n))
    p = lambda v: 1.0/(1.0 + np.exp(-v))
    return n, i1, j1, p(b1), i2, j2, p(b2), b2 - b1

def lowest_cusps(path, nmax):
    import pyarrow.csv as pc
    t = pc.read_csv(path, convert_options=pc.ConvertOptions(include_columns=["n", "i", "j", "pstar"]))
    n, i, j, p = (t[c].to_numpy() for c in ("n", "i", "j", "pstar"))
    m = n <= nmax; n, i, j, p = n[m], i[m], j[m], p[m]
    o = np.lexsort((p, n)); n, i, j, p = n[o], i[o], j[o], p[o]
    ns, first = np.unique(n, return_index=True)
    return {int(a): (int(i[k]), int(j[k]), float(p[k])) for a, k in zip(ns, first)}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="cusps/cusps_all.csv")
    ap.add_argument("--nmin", type=int, default=3)
    ap.add_argument("--nmax", type=int, default=5000)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--out", default="analysis")
    a = ap.parse_args()
    from multiprocessing import Pool
    _two_lowest(4, lnC_arr(4))                         # compile once before forking
    ns = list(range(a.nmin, a.nmax + 1))
    with Pool(a.workers) as pool:                      # large n first, for balance
        rows = sorted(pool.map(first_tie, sorted(ns, reverse=True), chunksize=4))
    fc = lowest_cusps(a.csv, a.nmax)
    out, same = [], []
    for n, i1, j1, p1, i2, j2, p2, gap in rows:
        ci, cj, cp = fc[n]
        s = (i1, j1) == (ci, cj)
        if s: same.append(n)
        out.append([n, i1, j1, j1 - i1, i1 + j1 - n, repr(float(p1)), f"{gap:.3e}",
                    ci, cj, cj - ci, ci + cj - n, repr(cp), int(s)])
    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, "first_tie_vs_cusp.csv")
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["n", "tie_i", "tie_j", "tie_width", "tie_band", "tie_pstar", "lnrho_gap_to_2nd",
                    "cusp_i", "cusp_j", "cusp_width", "cusp_band", "cusp_pstar", "same"])
        w.writerows(out)
    gaps = np.array([r[7] for r in rows])
    print(f"n = {a.nmin}..{a.nmax}: first tie point IS the first cusp at {len(same)} n: {same}")
    print(f"smallest ln-rho gap between first and second tie point: {gaps.min():.3e} "
          f"(n={rows[int(np.argmin(gaps))][0]})")
    print(path)

if __name__ == "__main__":
    main()
