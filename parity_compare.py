"""
parity_compare.py -- do odd and even n differ in their per-n cusp statistics?

    .venv/bin/python parity_compare.py [--nmin 100] [--out analysis/]

Reads public/per_n_summary.csv and analysis/prime_compare.csv (median and minimum of the cusp pair
mass f(i); run prime_compare.py first).  Nothing is recomputed from the cusp table.

SAWTOOTH TEST.  For each n, z(n) = s(n) - (s(n-1) + s(n+1))/2: n against the mean of its two
neighbours, which are the other parity.  The smooth trend cancels to second order, and a parity
offset delta (odd above even) gives z = +delta at odd n and -delta at even n, so
    delta = (mean z over odd n - mean z over even n)/2.
s is log(statistic) when the column is strictly positive (delta is then a relative offset,
reported in %), and the raw value otherwise (delta in the column's units, and as a fraction of the
column's local scatter).  p is two-sided, from 20,000 random reassignments of the parity labels.
"""
import argparse, csv, os
import numpy as np

COLS = ["n_cusps", "n_negF3", "n_double_cusps", "max_pstar", "min_pstar", "min_F3",
        "min_E_minus_Ehalf", "median_slope_jump", "min_slope_jump", "median_f", "min_f"]

def load():
    rows = {int(r["n"]): r for r in csv.DictReader(open("public/per_n_summary.csv"))}
    for r in csv.DictReader(open("analysis/prime_compare.csv")):
        rows[int(r["n"])].update(median_f=r["median_f"], min_f=r["min_f"])
    ns = np.array(sorted(rows))
    return ns, {c: np.array([float(rows[n][c]) for n in ns]) for c in COLS}

def sawtooth(ns, s, nmin, rng, reps=20000):
    z = s[1:-1] - (s[:-2] + s[2:])/2; m = ns[1:-1]
    k = m >= nmin; z, odd = z[k], (m[k] % 2 == 1)
    delta = (z[odd].mean() - z[~odd].mean())/2
    null = np.empty(reps)
    for t in range(reps):
        o = rng.permutation(odd)
        null[t] = (z[o].mean() - z[~o].mean())/2
    return delta, (np.abs(null) >= abs(delta)).mean(), z.std()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=100)
    ap.add_argument("--out", default="analysis")
    a = ap.parse_args()
    ns, data = load()
    assert np.array_equal(ns, np.arange(ns[0], ns[-1] + 1)), "gaps in n"
    rng = np.random.default_rng(1)
    print(f"n = {a.nmin}..{ns[-1]}: parity offset delta = odd minus even (sawtooth estimate)\n")
    print(f"{'statistic':<20} {'scale':<5} {'delta (odd-even)':>18} {'delta/scatter':>14} {'perm p':>8}")
    out = []
    for c in COLS:
        v = data[c]; lg = bool((v > 0).all())
        d, p, sd = sawtooth(ns, np.log(v) if lg else v, a.nmin, rng)
        shown = f"{100*(np.exp(d)-1):+.4f}%" if lg else f"{d:+.4g}"
        print(f"{c:<20} {'log' if lg else 'lin':<5} {shown:>18} {d/sd:>+14.3f} {p:>8.4f}")
        out.append([c, "log" if lg else "linear", f"{d:.6g}", f"{d/sd:.4g}", f"{p:.4g}"])
    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, "parity_compare.csv")
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["statistic", "scale", "delta_odd_minus_even", "delta_over_scatter", "perm_p"])
        w.writerows(out)
    print(f"\n{path}")

if __name__ == "__main__":
    main()
