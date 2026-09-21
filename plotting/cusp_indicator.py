"""
Plot the cusp indicator g = u(u+1), u = S_-/kappa, against tie-point index.

    .venv/bin/python plotting/cusp_indicator.py --n 1000 [--width 10000] [--height 6000]

h < 0 exactly at cusp points, so the zero line separates cusps from everything else:

    cusp  <=>  S_- < 0 < S_+        and    S_+ = S_- + kappa,  kappa = (j-i)f(i) > 0
          <=>  S_- * S_+ < 0        (kappa > 0 means S_+ > S_-, so only one ordering is possible)
          <=>  E'_- * E'_+ < 0      (the one-sided slopes; the derivative changes sign)

Writing u = S_-/kappa puts it in scale-free form: g = u(u+1), negative exactly on -1 < u < 0, with
a minimum of -1/4 at u = -1/2, where the zero sits dead centre in the slope jump.  But g reaches
1e304 at n=1000, so it cannot be plotted without throwing most of the range away.  Normalising by
magnitude instead keeps the sign and bounds the result:

    h = S_- S_+ / (S_-^2 + S_+^2)  in [-1/2, +1/2],   h < 0 <=> cusp

h = -1/2 is the perfectly symmetric cusp (S_- = -S_+, the zero dead centre in the jump); h -> +1/2
is as far from a cusp as a tie point gets.  The DEPTH below zero measures how robustly a point is a
cusp, not just whether it is one.

Everything is computed through ln_fi rather than kappa directly: kappa underflows to 0 in double
precision for deep tie points (12,792 of 249,001 at n=1000), which would send u to infinity.
"""
import argparse, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import cusps_data as cd

def indicator(d):
    """(sign(u), ln|u|) with u = S_-/kappa, via ln_fi so kappa may underflow.  Cusp <=> -1<u<0."""
    ln_kappa = np.log(d['j'].astype(float) - d['i'].astype(float)) + d['ln_fi']
    s = d['S_minus']
    with np.errstate(divide='ignore'):
        ln_abs_u = np.log(np.abs(s)) - ln_kappa
    return np.sign(s), ln_abs_u

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--data", default="data")
    ap.add_argument("--out", default="plots")
    ap.add_argument("--width", type=int, default=10000)
    ap.add_argument("--height", type=int, default=6000)
    ap.add_argument("--cap", type=float, default=1e4,
                    help="|u| is capped here for display; u reaches 1e300 and cannot be plotted raw")
    a = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    d = cd.load_ties(a.n, data=a.data)
    sign, ln_abs_u = indicator(d)
    cap = np.log(a.cap)
    u = sign*np.exp(np.minimum(ln_abs_u, cap))
    n_capped = int((ln_abs_u > cap).sum())
    x = d['rank_in_n']
    cusp = d['is_cusp']

    dpi = 100
    fig, ax = plt.subplots(figsize=(a.width/dpi, a.height/dpi), dpi=dpi)
    ax.axhline(0, color="0.35", lw=3, zorder=1)
    ax.axhspan(-1, 0, color="#d1495b", alpha=0.10, zorder=0)
    ax.scatter(x[~cusp], u[~cusp], s=1.5, c="#3b6ea8", alpha=0.35, linewidths=0,
               rasterized=True, zorder=2, label=f"tie points ({int((~cusp).sum()):,})")
    ax.scatter(x[cusp], u[cusp], s=90, c="#d1495b", linewidths=0,
               zorder=3, label=f"cusp points ({int(cusp.sum()):,})")
    ax.set_yscale("symlog", linthresh=1.0, linscale=3.0)
    ax.set_ylim(-a.cap*1.3, a.cap*1.3)
    ax.set_xlim(-0.005*len(x), 1.005*len(x))
    ax.set_xlabel("tie-point index, ordered by $p^*$ (0 = first tie point above 1/2)", fontsize=34, labelpad=25)
    ax.set_ylabel(r"$u=S_-/\kappa$   (position of zero within the slope jump)"
                  "\n" r"cusp $\Leftrightarrow\ -1<u<0$  (shaded)", fontsize=34, labelpad=25)
    ax.set_title(f"Cusp indicator across all {len(x):,} tie points of n={a.n}   "
                 f"(p* from {d['pstar'].min():.6f} to {d['pstar'].max():.6f})", fontsize=44, pad=35)
    ax.tick_params(labelsize=26, length=12, width=2)
    ax.grid(True, which="major", alpha=0.25, lw=1.2)
    leg = ax.legend(fontsize=30, markerscale=3, loc="upper left", framealpha=0.9)
    ax.text(0.995, 0.02, f"symlog: linear on |u|<1 (the cusp band), log outside; "
            f"|u| capped at {a.cap:g} for display ({n_capped:,} of {len(x):,} points, max |u| ~ 1e{ln_abs_u.max()/np.log(10):.0f})",
            transform=ax.transAxes, ha="right", fontsize=24, color="0.35")
    # second x axis in p*
    sec = ax.twiny(); sec.set_xlim(ax.get_xlim())
    ticks = np.linspace(0, len(x)-1, 11).astype(int)
    sec.set_xticks(ticks); sec.set_xticklabels([f"{d['pstar'][t]:.4f}" for t in ticks])
    sec.set_xlabel("$p^*$ at that index", fontsize=30, labelpad=20); sec.tick_params(labelsize=24)
    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, f"cusp_indicator_n{a.n:05d}.png")
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    print(f"{path}  ({os.path.getsize(path)/1e6:.1f} MB, {a.width}x{a.height})")
    print(f"  {len(x):,} tie points, {int(cusp.sum()):,} cusps; "
          f"u at cusps in [{u[cusp].min():.4f}, {u[cusp].max():.4f}]; "
          f"cusps span index {x[cusp].min():,}..{x[cusp].max():,} "
          f"(p* {d['pstar'][cusp].min():.5f}..{d['pstar'][cusp].max():.5f})")

if __name__ == "__main__":
    main()
