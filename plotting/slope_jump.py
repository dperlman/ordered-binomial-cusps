"""
The slope jump D = E'_+ - E'_- at every tie point, against tie-point index, coloured by pair mass.

    .venv/bin/python plotting/slope_jump.py --n 1000

D = (j-i) f(i) / (p* q*) > 0 at every tie point (RESEARCH_LOG.md section 2, fact 11: the kink at
a tie point is convex).  In this pipeline that positivity is built in -- S_+ is not stored but
derived as S_- + (j-i) f(i) -- so this plot cannot show a negative value; what it shows is the
MAGNITUDE, which spans ~300 orders, so y is log10 D computed in log space (kappa itself underflows
to 0 in double for the deepest ties).
"""
import argparse, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import cusps_data as cd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--data", default="data"); ap.add_argument("--out", default="plots")
    ap.add_argument("--width", type=int, default=10000); ap.add_argument("--height", type=int, default=6000)
    ap.add_argument("--cmin", type=float, default=-30)
    ap.add_argument("--linear", action="store_true", help="plot D itself instead of log10 D")
    ap.add_argument("--pmax", type=float, default=None, help="keep only tie points with p* <= pmax")
    ap.add_argument("--no-axis", action="store_true", help="leave out the p=1/2 axis row")
    ap.add_argument("--connect", action="store_true",
                    help="thin black line through consecutive tie points, drawn under the dots")
    a = ap.parse_args()
    import matplotlib; matplotlib.use("Agg")
    matplotlib.rcParams["agg.path.chunksize"] = 20000    # the --connect line can have >1M vertices
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.lines import Line2D

    d = cd.load_ties(a.n, data=a.data)
    keep = np.ones(len(d['pstar']), bool)
    if a.pmax is not None: keep &= d['pstar'] <= a.pmax
    if a.no_axis: keep &= ~d['is_axis']
    d = {k: (v[keep] if isinstance(v, np.ndarray) and v.shape[:1] == keep.shape else v) for k, v in d.items()}
    axis = d['is_axis']; cusp = d['is_cusp'] & ~axis; o = ~axis
    p = d['pstar']; w = d['j'].astype(float) - d['i'].astype(float)
    log10D = (np.log(w) + d['ln_fi'] - np.log(p*(1-p)))/np.log(10)    # exact, no underflow
    x = d['rank_in_n'].astype(float)
    col = np.clip(d['ln_fi']/np.log(10), a.cmin, 0.0); norm = Normalize(vmin=a.cmin, vmax=0.0)

    dpi = 100
    fig, ax = plt.subplots(figsize=(a.width/dpi, a.height/dpi), dpi=dpi)
    y = 10.0**log10D if a.linear else log10D
    ms = float(np.clip(1.2e6/o.sum(), 4.0, 60))          # bigger dots when there are few points
    if a.connect:
        order = np.argsort(x)
        ax.plot(x[order], y[order], color="black", lw=0.7, alpha=0.3, zorder=1, rasterized=True)
    sc = ax.scatter(x[o], y[o], c=col[o], s=ms, alpha=0.85, linewidths=0, cmap="viridis",
                    norm=norm, rasterized=True, zorder=2)
    ax.scatter(x[cusp], y[cusp], s=220, facecolors="none", edgecolors="#d1495b",
               linewidths=2.2, zorder=4)
    if axis.any():
        ax.scatter(x[axis], y[axis], s=1100, marker="*", c="#f0a202", edgecolors="#6b4500",
                   linewidths=3, zorder=5)
    ax.axhline(1.0 if a.linear else 0.0, color="0.3", lw=2.5, ls="--", zorder=1)
    if a.linear:
        ax.axhline(0.0, color="0.15", lw=3, zorder=1)
        top = max(y[o].max(), y[axis].max()) if axis.any() else y[o].max()
        ax.set_ylim(-0.03*top, 1.06*top)
    ax.set_xlim(-0.005*len(x), 1.005*len(x))
    ax.set_xlabel("tie-point index, ordered by $p^*$  (0 = the $p=1/2$ axis)", fontsize=36, labelpad=22)
    ax.set_ylabel((r"$D = E'_+ - E'_- = (j-i)\,f(i)/(p^*q^*)$   (linear)" if a.linear else
                   r"$\log_{10}\,D$,   $D = E'_+ - E'_- = (j-i)\,f(i)/(p^*q^*)$"), fontsize=36, labelpad=22)
    sec = ax.twiny(); sec.set_xlim(ax.get_xlim())
    ticks = np.linspace(0, len(x)-1, 11).astype(int)
    sec.set_xticks(ticks); sec.set_xticklabels([f"{p[t]:.4f}" for t in ticks])
    sec.set_xlabel("$p^*$ at that index", fontsize=30, labelpad=18); sec.tick_params(labelsize=24)
    rng = f" with $p^* \\leq {a.pmax:g}$" if a.pmax is not None else ""
    ax.set_title(f"Slope jump at the {int(o.sum()):,} tie points of n={a.n}{rng}, coloured by pair mass   "
                 f"(every value is positive: fact 11)", fontsize=44, pad=32)
    ax.tick_params(labelsize=26, length=12, width=2); ax.grid(True, alpha=0.22, lw=1.2)
    cb = fig.colorbar(sc, ax=ax, pad=0.012, fraction=0.025)
    cb.set_label(r"$\log_{10} f(i)$  (pair mass at the tie; $\leq$ %d clamped)" % a.cmin, fontsize=30, labelpad=18)
    cb.ax.tick_params(labelsize=24)
    handles = [Line2D([], [], marker="o", ls="", ms=14, mfc="none", mec="#d1495b", mew=2,
                      label=f"cusp ({int(cusp.sum())})"),
               *([Line2D([], [], marker="*", ls="", ms=26, mfc="#f0a202", mec="#6b4500", mew=2,
                        label=r"$p=1/2$ axis (all mirror pairs tie: $D = 2E'_+$)")] if axis.any() else []),
               Line2D([], [], ls="--", color="0.3", lw=2.5, label="$D = 1$")]
    ax.legend(handles=handles, fontsize=28, loc="upper right", framealpha=0.92)
    ax.text(0.995, 0.02, f"log10 D spans {log10D[o].min():.0f} .. {log10D[o].max():.2f};  "
            f"min D among cusps = 10^{log10D[cusp].min():.2f}",
            transform=ax.transAxes, ha="right", fontsize=24, color="0.35")
    os.makedirs(a.out, exist_ok=True)
    tag = ('_linear' if a.linear else '') + (f'_p{a.pmax:g}' if a.pmax is not None else '')
    path = os.path.join(a.out, f"slope_jump_n{a.n:05d}{tag}.png")
    fig.savefig(path, dpi=dpi, bbox_inches="tight"); plt.close(fig)
    print(f"{path}  ({os.path.getsize(path)/1e6:.1f} MB)   log10 D: min {log10D[o].min():.1f}, "
          f"max {log10D[o].max():.2f};  cusps span {log10D[cusp].min():.2f}..{log10D[cusp].max():.2f}")

if __name__ == "__main__":
    main()
