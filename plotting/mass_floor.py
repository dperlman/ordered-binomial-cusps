"""
Plot the pair mass f(i) at cusp points against n: the whole distribution, its per-n minimum,
and its per-n 1st percentile.

    .venv/bin/python plotting/mass_floor.py [--csv cusps/cusps_all.csv] [--data data]

A cusp needs |S_-| < kappa = (j-i) f(i): the kink has to be big enough to carry E' across zero, so
a tie point between two negligible masses cannot be a cusp however the rest of the ranking falls.
That puts a floor under f(i) at cusps.  This plot asks what the floor looks like as a function of n.

f(i) is recovered as kappa/(j-i) with kappa = S_+ - S_-.  That subtraction is catastrophic for deep
tie points -- it is why S_+ is not stored anywhere -- but at a CUSP kappa is comparable to |S_-| by
definition, so it is well conditioned here and good to ~10 digits.  The Parquet partitions carry
ln_fi directly, which is how the n >= 4000 points are obtained.

The per-n minimum is a near-miss statistic: it records how close the nearest tie point happened to
land to a smooth local maximum of E, so it scatters by orders of magnitude between neighbouring n.
The 1st percentile is drawn alongside, though with only 105-1063 cusps per n it is just the 1st to
11th smallest value and is itself noisy.

What is solid: the median and the maximum both scale as n^-1/2 with log-log correlation -1.00,
which is simply the 1/sqrt(n) scaling of the masses themselves.  What is not: the low tail follows
no law.  Fitting it gives n^-1.67 for the minimum (corr -0.55), n^-0.97 for the 1st percentile
(corr -0.27) and n^-3.32 for the 5th smallest -- three incompatible answers, because these are rare
events rather than a trend.  The running minimum is flat at 6.56e-8 from n=1075 through n=3000.
"""
import argparse, csv, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def from_csv(path):
    """{n: array of f(i) over that n's cusps} from the certified table."""
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            if r['pstar'] == '': continue
            n = int(r['n']); w = int(r['j']) - int(r['i'])
            out.setdefault(n, []).append((float(r['S_plus']) - float(r['S_minus']))/w)
    return {n: np.array(v) for n, v in out.items()}

def from_parquet(data, ns):
    """{n: array of f(i)} from the Parquet cusp partitions, excluding the p=1/2 axis row."""
    import cusps_data as cd
    out = {}
    for n in ns:
        try: d = cd.load_cusps(n, data=data)
        except FileNotFoundError: continue
        m = ~d['is_axis']
        out[n] = np.exp(d['ln_fi'][m])
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="cusps/cusps_all.csv")
    ap.add_argument("--data", default="data")
    ap.add_argument("--out", default="plots")
    ap.add_argument("--width", type=int, default=10000)
    ap.add_argument("--height", type=int, default=6000)
    ap.add_argument("--fit-from", type=int, default=300)
    a = ap.parse_args()
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    byn = from_csv(a.csv)
    ext = from_parquet(a.data, (4000, 5000, 6000, 7000, 8000))
    ns = np.array(sorted(byn))
    mn = np.array([byn[n].min() for n in ns])
    p1 = np.array([np.percentile(byn[n], 1) for n in ns])
    med = np.array([np.median(byn[n]) for n in ns])
    en = np.array(sorted(ext)); emn = np.array([ext[n].min() for n in en])
    ep1 = np.array([np.percentile(ext[n], 1) for n in en])

    dpi = 100
    fig, ax = plt.subplots(figsize=(a.width/dpi, a.height/dpi), dpi=dpi)
    # the full distribution, faintly
    allx = np.concatenate([np.full(len(byn[n]), n) for n in ns])
    ally = np.concatenate([byn[n] for n in ns])
    ax.scatter(allx, ally, s=1.0, c="#9bb7d4", alpha=0.18, linewidths=0, rasterized=True, zorder=1)
    # per-n statistics are separate values, not a continuous curve: draw them as points
    ax.plot(ns, med, lw=3.0, color="#2b5d8a", zorder=3)
    ax.scatter(ns, p1, s=9, c="#e08b2d", linewidths=0, alpha=0.85, rasterized=True, zorder=4)
    ax.scatter(ns, mn, s=9, c="#c1272d", linewidths=0, alpha=0.85, rasterized=True, zorder=5)
    # running lower envelope of the minimum
    env = np.minimum.accumulate(mn)
    ax.plot(ns, env, lw=3.5, color="#7a1216", ls="--", zorder=6)
    if len(en):
        ax.scatter(en, emn, s=700, marker="D", c="#c1272d", edgecolors="black", linewidths=2.5,
                   zorder=7)
        ax.scatter(en, ep1, s=700, marker="D", c="#e08b2d", edgecolors="black", linewidths=2.5,
                   zorder=7)
    # the BULK follows n^-1/2 cleanly (log-log corr -1.00); the low tail does not follow anything
    m = ns >= a.fit_from
    em, Am = np.polyfit(np.log(ns[m]), np.log(med[m]), 1)
    xf = np.array([a.fit_from, 8000.0])
    ax.plot(xf, np.exp(Am)*xf**em, lw=3, color="black", ls=":", zorder=8,
            )
    e, A = np.polyfit(np.log(ns[m]), np.log(mn[m]), 1)
    rmin = np.corrcoef(np.log(ns[m]), np.log(mn[m]))[0,1]
    lo = np.array([byn[n].min() for n in ns]).min()
    k = ns[np.argmin(mn)]
    ax.annotate(f"smallest pair mass at any cusp, n $\\leq$ 3000:\n{lo:.2e}  (n={k})",
                xy=(k, lo), xytext=(k*1.6, lo*0.22), fontsize=28,
                arrowprops=dict(arrowstyle="->", lw=2.5, color="0.3"))
    ax.text(0.615, 0.975,
            "The BULK scales as $n^{-1/2}$ (median and max both, log-log corr $-1.00$): the masses\n"
            "themselves shrink that way.  The LOW TAIL does not follow a law -- fitting it gives\n"
            f"$n^{{{e:.2f}}}$ for the minimum (corr {rmin:.2f}), $n^{{-0.97}}$ for the 1st percentile (corr $-0.27$),\n"
            "$n^{-3.32}$ for the 5th smallest.  Those are near-miss events, not a trend.\n"
            "The 1st percentile is only the 1st-11th smallest value here (105-1063 cusps per $n$).",
            transform=ax.transAxes, va="top", fontsize=26, color="0.25",
            bbox=dict(boxstyle="round,pad=0.6", fc="white", ec="0.7", alpha=0.93))
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("$n$", fontsize=38, labelpad=22)
    ax.set_ylabel(r"pair mass at the tie,  $f(i)=f(j)$", fontsize=38, labelpad=22)
    ax.set_title("How small can the tying masses be at a cusp?   "
                 r"(a cusp needs $|S_-| < \kappa = (j-i)\,f(i)$)", fontsize=44, pad=34)
    ax.tick_params(labelsize=28, length=12, width=2, which="major")
    ax.tick_params(length=6, width=1.2, which="minor")
    ax.grid(True, which="major", alpha=0.28, lw=1.2)
    ax.grid(True, which="minor", alpha=0.12, lw=0.8)
    from matplotlib.lines import Line2D
    dot = lambda c, l, ms=13: Line2D([], [], marker="o", ls="", ms=ms, mfc=c, mec=c, label=l)
    handles = [dot("#9bb7d4", f"every cusp ({len(ally):,})", 11),
               Line2D([], [], lw=3.0, color="#2b5d8a", label="median over each n"),
               dot("#e08b2d", "1st percentile over each n"),
               dot("#c1272d", "minimum over each n"),
               Line2D([], [], lw=3.5, ls="--", color="#7a1216", label="running minimum (lower envelope)"),
               Line2D([], [], marker="D", ls="", ms=15, mfc="#c1272d", mec="black", mew=1.8,
                      label="minimum, n = 4000…8000 (Parquet)"),
               Line2D([], [], marker="D", ls="", ms=15, mfc="#e08b2d", mec="black", mew=1.8,
                      label="1st percentile, n = 4000…8000"),
               Line2D([], [], lw=3, ls=":", color="black",
                      label=f"fit to the median:  $f \\sim {np.exp(Am):.2f}\\,n^{{{em:.2f}}}$   "
                            f"(corr {np.corrcoef(np.log(ns[m]), np.log(med[m]))[0,1]:.2f})")]
    ax.legend(handles=handles, fontsize=27, loc="lower left", framealpha=0.93)
    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, "cusp_mass_floor.png")
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    print(f"{path}  ({os.path.getsize(path)/1e6:.1f} MB)")
    print(f"  n = {ns.min()}..{ns.max()} from the CSV, plus {list(en)} from Parquet")
    print(f"  overall minimum {lo:.3e} at n={k}")
    print(f"  median ~ {np.exp(Am):.3f} n^{em:.3f} (corr {np.corrcoef(np.log(ns[m]), np.log(med[m]))[0,1]:.3f});  "
          f"minimum ~ n^{e:.2f} (corr {rmin:.2f}, unreliable)")
    print(f"  1st percentile at n=1000/2000/3000: "
          f"{p1[ns==1000][0]:.2e} / {p1[ns==2000][0]:.2e} / {p1[ns==3000][0]:.2e}")
    if len(en): print(f"  minimum at n=4000..8000: {', '.join(f'{v:.2e}' for v in emn)}")

if __name__ == "__main__":
    main()
