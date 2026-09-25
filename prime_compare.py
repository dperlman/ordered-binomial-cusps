"""
prime_compare.py -- do prime n (or n with n+1 prime) differ in their cusp statistics?

    .venv/bin/python prime_compare.py [--csv cusps/cusps_all.csv] [--nmin 100] [--out analysis/]

Per n: the cusp count, and the median and minimum of the pair mass f(i) = (S_+ - S_-)/(j-i).

PARITY IS THE CONFOUNDER.  Every prime above 2 is odd, and n+1 prime (n >= 2) forces n even, while
parity alone splits some per-n statistics (plotting/lowest_cusp.py: E - E(1/2) of the lowest cusp
sits on two separate curves for even and odd n).  So the comparisons are made WITHIN a parity:
    n prime      vs  odd composite n
    n+1 prime    vs  even n with n+1 composite
Each statistic is detrended locally: residual(n) = log s(n) - median of log s over the SAME-parity
n in [n-W, n+W], n itself excluded.  The group difference is the mean residual of one group minus
the other, with a two-sided permutation p-value (labels shuffled within the parity class).
A parity row (odd vs even, detrended over all n) is printed first as a control (parity_compare.py
tests parity over every per-n column).

Writes analysis/prime_compare.csv (one row per n: stats, residuals, flags).  numpy + pyarrow.
"""
import argparse, csv, os
import numpy as np

def per_n(path):
    import pyarrow.csv as pc
    t = pc.read_csv(path, convert_options=pc.ConvertOptions(
        include_columns=["n", "i", "j", "S_minus", "S_plus"]))
    n = t["n"].to_numpy()
    f = (t["S_plus"].to_numpy() - t["S_minus"].to_numpy())/(t["j"].to_numpy() - t["i"].to_numpy())
    o = np.argsort(n, kind="stable"); n, f = n[o], f[o]
    ns, start, cnt = np.unique(n, return_index=True, return_counts=True)
    g = np.split(f, start[1:])
    return ns, cnt, np.array([np.median(x) for x in g]), np.array([x.min() for x in g])

def sieve(top):
    s = np.ones(top + 1, bool); s[:2] = False
    for q in range(2, int(top**0.5) + 1):
        if s[q]: s[q*q::q] = False
    return s

def residual(ns, y, W, same_parity=True):
    """log y(n) minus the median of log y over neighbours within W (same parity), self excluded."""
    ly = np.log(y); r = np.empty(len(ns))
    for a, n in enumerate(ns):
        m = (np.abs(ns - n) <= W) & (ns != n)
        if same_parity: m &= (ns % 2 == n % 2)
        r[a] = ly[a] - np.median(ly[m])
    return r

def compare(r, grp, rng, reps=20000):
    """(n1, n0, mean diff, perm p) for residuals r split by boolean grp."""
    d = r[grp].mean() - r[~grp].mean()
    k = grp.sum(); tot = r.sum(); null = np.empty(reps)
    for t in range(reps):
        s = r[rng.permutation(len(r))[:k]].sum()
        null[t] = s/k - (tot - s)/(len(r) - k)
    return k, len(r) - k, d, (np.abs(null) >= abs(d)).mean()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="cusps/cusps_all.csv")
    ap.add_argument("--nmin", type=int, default=100)
    ap.add_argument("--window", type=int, default=100)
    ap.add_argument("--out", default="analysis")
    a = ap.parse_args()
    ns, cnt, med, mn = per_n(a.csv)
    P = sieve(int(ns.max()) + 1)
    stats = {"count": cnt.astype(float), "median f": med, "min f": mn}
    res = {k: residual(ns, v, a.window) for k, v in stats.items()}
    keep = ns >= a.nmin
    odd, even = keep & (ns % 2 == 1), keep & (ns % 2 == 0)
    rng = np.random.default_rng(1)
    print(f"n = {a.nmin}..{ns.max()}, local window +-{a.window} (same parity), "
          f"differences in mean log-residual (x100 = ~percent)\n")
    print(f"{'statistic':<10} {'comparison':<32} {'n_yes':>6} {'n_no':>6} {'diff':>9} {'perm p':>8}")
    for k, v in stats.items():
        rp = residual(ns, v, a.window, same_parity=False)
        n1, n0, d, p = compare(rp[keep], (ns % 2 == 1)[keep], rng)
        print(f"{k:<10} {'odd vs even (control)':<32} {n1:>6} {n0:>6} {d:>+9.4f} {p:>8.4f}")
        for lab, cls, flag in (("n prime vs odd composite", odd, P[ns]),
                               ("n+1 prime vs even, n+1 comp.", even, P[ns + 1])):
            n1, n0, d, p = compare(res[k][cls], flag[cls], rng)
            print(f"{'':<10} {lab:<32} {n1:>6} {n0:>6} {d:>+9.4f} {p:>8.4f}")
    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, "prime_compare.csv")
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["n", "n_prime", "n1_prime", "n_cusps", "median_f", "min_f",
                    "resid_count", "resid_median_f", "resid_min_f"])
        for a_, n in enumerate(ns):
            w.writerow([n, int(P[n]), int(P[n+1]), cnt[a_],
                        repr(float(med[a_])), repr(float(mn[a_])), f"{res['count'][a_]:.6g}", f"{res['median f'][a_]:.6g}",
                        f"{res['min f'][a_]:.6g}"])
    print(f"\n{path}")

if __name__ == "__main__":
    main()
