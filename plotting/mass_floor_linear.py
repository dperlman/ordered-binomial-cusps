"""plotting/mass_floor_linear.py -- the pair mass f(i) at cusp points against n, LINEAR n axis,
one pixel column per n (plotting/_style.py, NGrid).

    .venv/bin/python plotting/mass_floor_linear.py [--csv cusps/cusps_all.csv] [--out plots]

The linear-axis companion of mass_floor.py (log-log, same quantity; its docstring has the
reasoning).  Every n of the complete table -- n=3..5000, no Parquet extension -- owns exactly one
column.  The cloud is every cusp binned into its pixel and shaded by count; the per-n minimum and
median are 1 px wide marks in that same column; the running minimum and the n^-1/2 fit to the
median are drawn over them.  (A prime-n strip and the 1st percentile were tried and dropped
2026-09-23: neither added anything readable.)  n starts at 3 because n=2 has no cusp: its only
tie point, (1,2) at p*=2/3, has left slope exactly 0.

f(i) = (S_+ - S_-)/(j-i): well conditioned at a cusp, where kappa is comparable to |S_-|.
"""
import argparse, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _style as st

C_MED, C_MIN, C_ENV = "#2a78d6", "#e34948", "#7a1216"   # from a validated set
CLOUD_LO, CLOUD_HI = "#dde1e6", "#3d434b"

def load(path):
    """(n, f) for every cusp, sorted by n."""
    import pyarrow.csv as pc
    t = pc.read_csv(path, convert_options=pc.ConvertOptions(
        include_columns=["n", "i", "j", "S_minus", "S_plus"]))
    n = t["n"].to_numpy(); w = t["j"].to_numpy() - t["i"].to_numpy()
    f = (t["S_plus"].to_numpy() - t["S_minus"].to_numpy())/w
    o = np.argsort(n, kind="stable")
    return n[o], f[o]

def per_n(n, f):
    ns, start = np.unique(n, return_index=True)
    groups = np.split(f, start[1:])
    return (ns, np.array([g.min() for g in groups]), np.array([np.median(g) for g in groups]))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="cusps/cusps_all.csv")
    ap.add_argument("--out", default="plots")
    ap.add_argument("--height", type=int, default=st.DEFAULT_H)
    ap.add_argument("--fit-from", type=int, default=300)
    a = ap.parse_args()

    n, f = load(a.csv)
    ns, mn, med = per_n(n, f)
    if not np.array_equal(ns, np.arange(ns[0], ns[-1] + 1)):
        sys.exit("the table has gaps in n; this plot is for a complete range")
    env = np.minimum.accumulate(mn)
    m = ns >= a.fit_from
    lx = np.log(ns[m])
    em, Am = np.polyfit(lx, np.log(med[m]), 1); rmed = np.corrcoef(lx, np.log(med[m]))[0, 1]
    e0, _ = np.polyfit(lx, np.log(mn[m]), 1);  rmin = np.corrcoef(lx, np.log(mn[m]))[0, 1]
    kmin = ns[np.argmin(mn)]; lo = mn.min()
    flat = ns[env == lo]                        # where the running minimum sits at its floor

    g = st.NGrid(ns[0], ns[-1], height=a.height)
    ax, back = g.ax, g.back[0]
    ax.set_yscale("log")
    ax.set_ylim(10**np.floor(np.log10(f.min()) - 0.3), 10**np.ceil(np.log10(f.max()) + 0.1))

    cmax = g.density(n, f, lo=CLOUD_LO, hi=CLOUD_HI)
    g.points(ns, med, C_MED, h=5)
    g.points(ns, mn, C_MIN, h=5)
    ax.plot(ns, env, drawstyle="steps-mid", lw=3.5, color=C_ENV, ls="--")
    xf = np.array([a.fit_from, ns[-1]], float)
    ax.plot(xf, np.exp(Am)*xf**em, lw=3, color="black", ls=":")

    ax.annotate(f"smallest pair mass at any cusp, n $\\leq$ {ns[-1]}:  {lo:.2e}  (n = {kmin})",
                xy=(kmin, lo), xytext=(kmin + 0.04*len(ns), lo*12), fontsize=28,
                arrowprops=dict(arrowstyle="->", lw=2.5, color="0.3"))
    ax.text(0.985, 0.975,
            f"The BULK scales as $n^{{{em:.2f}}}$ (median, fit from n = {a.fit_from}, "
            f"log-log corr {rmed:.2f}):\nthe masses themselves shrink that way.  The LOW TAIL "
            f"follows no law:\nfitting the minimum gives $n^{{{e0:.2f}}}$ with corr only {rmin:.2f}.\n"
            f"The running minimum stays at {lo:.2e} from n = {flat[0]} to {flat[-1]}.",
            transform=ax.transAxes, ha="right", va="top", fontsize=26, color="0.2",
            bbox=dict(boxstyle="round,pad=0.6", fc="white", ec="0.7", alpha=0.93))

    back.set_xlabel("$n$", fontsize=38, labelpad=18)
    back.set_ylabel(r"pair mass at the tie,  $f(i)=f(j)$", fontsize=38, labelpad=22)
    back.set_title("How small can the tying masses be at a cusp?   "
                   r"(a cusp needs $|S_-| < \kappa = (j-i)\,f(i)$)", fontsize=44, pad=30)
    back.tick_params(labelsize=28, length=12, width=2, which="major")
    back.tick_params(length=6, width=1.2, which="minor")
    back.grid(True, which="major", alpha=0.35, lw=1.2)
    back.grid(True, which="minor", axis="y", alpha=0.15, lw=0.8)
    back.set_axisbelow(True)

    mark = lambda c, l: Line2D([], [], marker="s", ls="", ms=13, mfc=c, mec=c, label=l)
    handles = [Patch(fc=CLOUD_HI, ec="none",
                     label=f"every cusp ({len(f):,}); darker = more per pixel (max {cmax})"),
               mark(C_MED, "median over each n"),
               mark(C_MIN, "minimum over each n"),
               Line2D([], [], lw=3.5, ls="--", color=C_ENV, label="running minimum"),
               Line2D([], [], lw=3, ls=":", color="black",
                      label=f"fit to the median:  $f \\sim {np.exp(Am):.2f}\\,n^{{{em:.2f}}}$")]
    ax.legend(handles=handles, fontsize=27, loc="lower left", framealpha=0.93)

    os.makedirs(a.out, exist_ok=True)
    path = g.save(os.path.join(a.out, "cusp_mass_floor_linear.png"))
    print(f"{path}  {g.W}x{g.H} px, data area {g.pw}x{g.ph} px, {g.k} px per n  "
          f"({os.path.getsize(path)/1e6:.1f} MB)")
    print(f"  n = {ns[0]}..{ns[-1]}, {len(f):,} cusps, max {cmax} per pixel")
    print(f"  minimum {lo:.3e} at n={kmin}; running minimum flat over n={flat[0]}..{flat[-1]}")
    print(f"  median ~ {np.exp(Am):.3f} n^{em:.3f} (corr {rmed:.3f}); "
          f"min ~ n^{e0:.2f} (corr {rmin:.2f})")

if __name__ == "__main__":
    main()
