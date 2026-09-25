"""
star_analyze.py -- analyses of the (★) switch-point check written by star_check.py.

    .venv/bin/python star_analyze.py --out analysis/ [--cusps cusps/cusps_all.csv]

Reads <out>/star_per_n.csv, <out>/star_switch_points/ (first 20 rows of each file) and the band-1..3
cusps of cusps_all.csv.  Writes <out>/star_bands.csv (one row per n: certificate margin and true cusp
margin for bands 1-3) and prints:
  1. sign of min_M over all n, and where the minimum over m sits; whether M(n,m) is strictly
     increasing in m (every switch point, all files); M / (2 sqrt(2/pi) x_m^2) and max |x_m - m/2|;
  2. retained fraction M(n,b) / T(n,b), T = (E(cusp) - E(1/2)) n^1.5 for the lowest cusp of band b
     (i+j = n+b), by parity and by n-range;
  3. the asymptotic of item 4c: R(n,m) = M(n,m) - 2 sqrt(2/pi) (m/2)^2 fitted as L + b/n + c/n^2
     over n in [500, 5000], separately for even and odd n, m = 1..20;
  4. runtime per n (star_timings.csv) and a power-law fit for extrapolation.
"""
import argparse, csv, glob, math, os
from multiprocessing import Pool
import numpy as np

import binom_core as core

def load_per_n(out):
    with open(os.path.join(out, "star_per_n.csv")) as fh:
        return list(csv.DictReader(fh))

def load_M20(out, ns, mmax=20):
    """M(n, m) for m = 1..mmax as an array indexed [n, m-1] (nan where m > n-1)."""
    A = np.full((max(ns) + 1, mmax), np.nan)
    for n in ns:
        with open(os.path.join(out, "star_switch_points", f"n{n:05d}.csv")) as fh:
            r = csv.reader(fh); next(r)
            for k, row in zip(range(mmax), r):
                A[n, int(row[1]) - 1] = float(row[7])
    return A

def band_cusps(path, nmax):
    """{(n, b): (E_min, i, j, count)} for b = 1..3, lowest-E cusp of each band."""
    import pyarrow.csv as pc
    t = pc.read_csv(path, convert_options=pc.ConvertOptions(include_columns=["n", "i", "j", "E"]))
    n, i, j, E = (t[c].to_numpy() for c in ("n", "i", "j", "E"))
    b = i + j - n
    sel = (b <= 3) & (n <= nmax)
    d = {}
    for nn, ii, jj, bb, ee in zip(n[sel], i[sel], j[sel], b[sel], E[sel]):
        key = (int(nn), int(bb))
        if key not in d or ee < d[key][0]: d[key] = (float(ee), int(ii), int(jj), d.get(key, (0, 0, 0, 0))[3] + 1)
        else: d[key] = d[key][:3] + (d[key][3] + 1,)
    return d

def whole_file(path):
    """(n, M strictly increasing in m?, min over m of M/(2 sqrt(2/pi) x^2) and its m, max |x_m - m/2|)."""
    import pyarrow.csv as pc
    t = pc.read_csv(path, convert_options=pc.ConvertOptions(include_columns=["n", "m", "x_m", "M"]))
    n = int(t["n"][0].as_py()); m = t["m"].to_numpy(); x = t["x_m"].to_numpy(); M = t["M"].to_numpy()
    r = M/(2*math.sqrt(2/math.pi)*x**2)
    return n, bool(np.all(np.diff(M) > 0)), float(r.min()), int(m[np.argmin(r)]), float(np.max(np.abs(x - m/2)))

def fit_LR(ns, R):
    """Least squares R = L + b/n + c/n^2; returns (L, b, c, rms)."""
    X = np.column_stack([np.ones_like(ns, dtype=float), 1.0/ns, 1.0/ns**2])
    coef, *_ = np.linalg.lstsq(X, R, rcond=None)
    return (*coef, float(np.sqrt(np.mean((X @ coef - R)**2))))

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="analysis")
    ap.add_argument("--cusps", default="cusps/cusps_all.csv")
    ap.add_argument("--fit-nmin", type=int, default=500)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    rows = load_per_n(a.out)
    ns = np.array([int(r["n"]) for r in rows]); nmax = int(ns.max())
    minM = np.array([float(r["min_M"]) for r in rows])
    arg = np.array([int(r["argmin_m"]) for r in rows])

    print("== 1. sign of the margin")
    print(f"  n = {ns.min()}..{nmax}: min_M > 0 for {np.sum(minM > 0)} of {len(ns)};  "
          f"smallest min_M = {minM.min():.6f} at n = {ns[np.argmin(minM)]}")
    for par, name in ((0, "even"), (1, "odd")):
        s = (ns % 2 == par) & (ns >= 100)
        print(f"  {name} n>=100: min_M in [{minM[s].min():.6f}, {minM[s].max():.6f}]")
    vals, cnt = np.unique(arg, return_counts=True)
    print(f"  argmin_m distribution: {dict(zip(vals.tolist(), cnt.tolist()))}")
    print(f"  n with argmin_m != 1: {ns[arg != 1].tolist()}")

    with Pool(a.workers) as pool:
        W = pool.map(whole_file, sorted(glob.glob(os.path.join(a.out, "star_switch_points", "n*.csv"))))
    notinc = [w[0] for w in W if not w[1]]
    print(f"  n where M(n,m) is NOT strictly increasing in m: {notinc[:20]} ({len(notinc)} of {len(W)})")
    for lo in (3, 100, 1000):
        r, n_, m_ = min((w[2], w[0], w[3]) for w in W if w[0] >= lo)
        print(f"  n>={lo}: min over n,m of M/(2 sqrt(2/pi) x_m^2) = {r:.4f} (n={n_}, m={m_})")
    print(f"  max |x_m - m/2| over every switch point: {max(w[4] for w in W):.4f} "
          f"(n={max((w[4], w[0]) for w in W)[1]})")

    A = load_M20(a.out, ns.tolist())
    C = 2*math.sqrt(2/math.pi)

    print("\n== 2. retained fraction against the true band-b cusp margin")
    bc = band_cusps(a.cusps, nmax)
    with open(os.path.join(a.out, "star_bands.csv"), "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["n", "b", "M", "T", "ratio", "cusp_i", "cusp_j", "n_cusps_in_band"])
        R = {}
        for n in ns.tolist():
            Eh = core.E_half(n)
            for b in (1, 2, 3):
                if (n, b) not in bc or b > n - 1: continue
                E, i, j, k = bc[(n, b)]
                T = (E - Eh)*n**1.5; M = float(A[n, b - 1])
                R[(n, b)] = M/T
                w.writerow([n, b, repr(M), repr(T), repr(M/T), i, j, k])
    bands_missing = {b: [n for n in ns.tolist() if (n, b) not in R] for b in (1, 2, 3)}
    print(f"  n with no cusp in band b: " + "; ".join(f"b={b}: {v[:12]}{'...' if len(v) > 12 else ''} ({len(v)})"
                                                   for b, v in bands_missing.items()))
    edges = [100, 500, 1000, 2000, 3000, 4000, 5001]
    for b in (1, 2, 3):
        for par, name in ((0, "even"), (1, "odd")):
            parts = []
            for lo, hi in zip(edges[:-1], edges[1:]):
                v = np.array([R[(n, b)] for n in range(lo, min(hi, nmax + 1)) if n % 2 == par and (n, b) in R])
                if len(v): parts.append(f"{lo}-{hi-1}: {np.median(v):.4f} [{v.min():.4f},{v.max():.4f}]")
            print(f"  band {b} {name:4s}: " + " | ".join(parts))

    print(f"\n== 3. R(n,m) = M(n,m) - 2 sqrt(2/pi) (m/2)^2 ~ L + b/n + c/n^2, n in [{a.fit_nmin}, {nmax}]")
    print(f"  predicted L: -+ 1/(4 sqrt(2 pi)) = -+{1/(4*math.sqrt(2*math.pi)):.5f} for odd m (minus n even), 0 for even m")
    print(f"  {'m':>3} {'L_even':>10} {'b_even':>9} {'c_even':>9} {'rms':>8} | {'L_odd':>10} {'b_odd':>9} {'c_odd':>9} {'rms':>8}"
          f" | {'(Le+Lo)/2':>10} {'(Lo-Le)/2':>10}")
    for m in range(1, 21):
        out = []
        for par in (0, 1):
            s = (ns >= a.fit_nmin) & (ns % 2 == par)
            nn = ns[s].astype(float); r = A[ns[s], m - 1] - C*(m/2)**2
            out.append(fit_LR(nn, r))
        (Le, be, ce, re), (Lo, bo, co, ro) = out
        print(f"  {m:3d} {Le:10.6f} {be:9.4f} {ce:9.1f} {re:8.1e} | {Lo:10.6f} {bo:9.4f} {co:9.1f} {ro:8.1e}"
              f" | {(Le+Lo)/2:10.6f} {(Lo-Le)/2:10.6f}")

    print("\n== 4. runtime")
    with open(os.path.join(a.out, "star_timings.csv")) as fh:
        T = np.array([(int(r["n"]), float(r["seconds"])) for r in csv.DictReader(fh)])
    s = T[:, 0] >= 1000
    k, lc = np.polyfit(np.log(T[s, 0]), np.log(T[s, 1]), 1)
    print(f"  per-n worker time ~ n^{k:.2f} over n>=1000;  t(1000) = {math.exp(lc)*1000**k:.3f}s, "
          f"t({nmax}) = {math.exp(lc)*nmax**k:.3f}s (measured {T[T[:,0]==nmax,1][0]:.3f}s)")
    print(f"  total worker time {T[:,1].sum():.0f}s over n=3..{nmax}")
    tN = T[T[:, 0] == nmax, 1][0]
    for N in (10000, 20000, 50000):
        extra = sum(tN*(n/nmax)**k for n in range(nmax + 1, N + 1))
        print(f"  extending to n={N}: ~{extra:.0f} CPU-s more, ~{extra/8/60:.1f} min on 8 workers "
              f"(t(N) ~ {tN*(N/nmax)**k:.1f}s; anchored on the measured t({nmax}))")
