"""
Plot the one-sided slopes E'_- and E'_+ at every tie point of one n, coloured by the pair mass.

    .venv/bin/python plotting/slope_plot.py --n 1000 --x p          # x = p*, linear
    .venv/bin/python plotting/slope_plot.py --n 1000 --x index      # x = tie-point index
    .venv/bin/python plotting/slope_plot.py --n 1000 --x both

We only know E' at the tie points, so each p* carries two dots: E'_- just before the swap and
E'_+ just after.  The vertical gap between them is the kink D = kappa/(p*q*); for most tie points
it is far below one pixel and the two dots coincide.  No lines are drawn -- we have no information
about E' between tie points beyond concavity.

Colour is log10 of the pair mass f(i) = f(j) at the tie.  It sets the kink (kappa = (j-i) f(i),
and f spans ~300 orders while j-i spans three, so corr(ln f, ln kappa) = 1.000), i.e. it colours
by how much the tie point can bend E.  Width j-i would give almost the same picture, distorted by
the p-dependence ln f ~ c - (j-i)^2/(8 n p q): "nearby" is measured in standard deviations, which
shrink toward p=1.
"""
import argparse, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import cusps_data as cd

def render(d, xmode, a):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    ax_row = d['is_axis']; cusp = d['is_cusp'] & ~ax_row; ordinary = ~ax_row
    x = d['pstar'] if xmode == "p" else d['rank_in_n'].astype(float)
    lo, hi = a.ymin, a.ymax
    sl = np.clip(d['slope_left'], lo, hi); sr = np.clip(d['slope_right'], lo, hi)
    n_pin = int(((d['slope_left'] > hi) | (d['slope_right'] > hi)).sum())
    col = np.clip(d['ln_fi']/np.log(10), a.cmin, 0.0)
    norm = Normalize(vmin=a.cmin, vmax=0.0)

    dpi = 100
    fig, ax = plt.subplots(figsize=(a.width/dpi, a.height/dpi), dpi=dpi)
    ax.axhline(0, color="0.25", lw=3, zorder=1)
    kw = dict(s=2.2, alpha=0.45, linewidths=0, cmap=a.cmap, norm=norm, rasterized=True)
    ax.scatter(x[ordinary], sl[ordinary], c=col[ordinary], zorder=2, **kw)
    sc = ax.scatter(x[ordinary], sr[ordinary], c=col[ordinary], zorder=2, **kw)
    # cusps: hollow rings on both dots of the straddling pair
    for y in (sl, sr):
        ax.scatter(x[cusp], y[cusp], s=140, facecolors="none", edgecolors="#d1495b",
                   linewidths=1.8, zorder=4)
    # the p=1/2 axis: two antisymmetric dots
    for y in (sl, sr):
        ax.scatter(x[ax_row], y[ax_row], s=1100, marker="*", c="#f0a202", edgecolors="#6b4500",
                   linewidths=3, zorder=5)
    ax.set_ylim(lo - 0.02*(hi-lo), hi + 0.02*(hi-lo))
    if xmode == "p":
        ax.set_xlim(0.5 - 0.005, 1.0 + 0.005)
        ax.set_xlabel("$p^*$", fontsize=36, labelpad=22)
    else:
        ax.set_xlim(-0.005*len(x), 1.005*len(x))
        ax.set_xlabel("tie-point index, ordered by $p^*$  (0 = the $p=1/2$ axis)", fontsize=36, labelpad=22)
        sec = ax.twiny(); sec.set_xlim(ax.get_xlim())
        ticks = np.linspace(0, len(x)-1, 11).astype(int)
        sec.set_xticks(ticks); sec.set_xticklabels([f"{d['pstar'][t]:.4f}" for t in ticks])
        sec.set_xlabel("$p^*$ at that index", fontsize=30, labelpad=18); sec.tick_params(labelsize=24)
    ax.set_ylabel(r"$E'(p)$ at tie points:  $E'_-$ and $E'_+$, one dot each", fontsize=36, labelpad=22)
    n = int(d['n'][0])
    ax.set_title(f"One-sided slopes of $E(n,p)$ at all {int(ordinary.sum()):,} tie points of n={n}, "
                 f"coloured by pair mass", fontsize=44, pad=32)
    ax.tick_params(labelsize=26, length=12, width=2)
    ax.grid(True, alpha=0.22, lw=1.2)
    cb = fig.colorbar(sc, ax=ax, pad=0.012, fraction=0.025)
    cb.set_label(r"$\log_{10} f(i)$  (pair mass at the tie; $\leq$ %d clamped)" % a.cmin, fontsize=30, labelpad=18)
    cb.ax.tick_params(labelsize=24)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker="o", ls="", ms=14, mfc="none", mec="#d1495b", mew=2,
                      label=f"cusp: the pair straddles zero ({int(cusp.sum())})"),
               Line2D([], [], marker="*", ls="", ms=26, mfc="#f0a202", mec="#6b4500", mew=2,
                      label=r"$p=1/2$ axis: $E'_\pm = \pm%.2f$" % d['slope_right'][ax_row][0])]
    ax.legend(handles=handles, fontsize=28, loc="upper left", framealpha=0.92)
    ax.text(0.995, 0.02, f"y clipped to [{lo:g}, {hi:g}]: {n_pin:,} dots pinned at the top edge "
            f"(max E' = {d['slope_right'].max():.0f}, near p=1 where p*q*"r"$\to$0)",
            transform=ax.transAxes, ha="right", fontsize=24, color="0.35")
    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, f"slope_n{n:05d}_{xmode}.png")
    fig.savefig(path, dpi=dpi, bbox_inches="tight"); plt.close(fig)
    print(f"{path}  ({os.path.getsize(path)/1e6:.1f} MB)  pinned {n_pin:,}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--x", choices=["p", "index", "both"], default="both")
    ap.add_argument("--data", default="data"); ap.add_argument("--out", default="plots")
    ap.add_argument("--width", type=int, default=10000); ap.add_argument("--height", type=int, default=6000)
    ap.add_argument("--ymin", type=float, default=-30); ap.add_argument("--ymax", type=float, default=120)
    ap.add_argument("--cmin", type=float, default=-30, help="log10 f(i) floor for the colour scale")
    ap.add_argument("--cmap", default="viridis")
    a = ap.parse_args()
    d = cd.load_ties(a.n, data=a.data)
    for m in (["p", "index"] if a.x == "both" else [a.x]): render(d, m, a)

if __name__ == "__main__":
    main()
