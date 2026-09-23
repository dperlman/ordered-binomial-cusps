"""plotting/lowest_cusp.py -- the lowest (minimum-E) cusp of each n, versus n.

Three figures, all with n on the x axis, from public/per_n_summary.csv (columns argmin_i/argmin_j
give the minimising cusp, min_E_minus_Ehalf its height above E(1/2)):

  lowest_cusp_E.png      E(p*) - E(1/2) at that cusp        (log-log, plus the n^-3/2 scaling)
  lowest_cusp_p.png      p* of that cusp                    (log-log, plus the 1/(2n) scaling)
  lowest_cusp_width.png  j - i of that cusp                 (the colouring quantity)

COLOURING.  The obvious choice was the BAND, i+j-n.  It carries no information: the lowest cusp is
in band 1 for all 4998 values of n.  Band 1 is a whole family though -- i+j = n+1 with any width
w = j-i -- and given band 1, w has the same parity as n+1, so w's parity is locked to n's.  That
parity is what actually splits the data: E-E(1/2) times n^1.5 sits at 0.4834 for even n and 0.6828
for odd n, two curves with no overlap at all.  So the points are coloured by n parity, and the
third figure plots the width itself.

Run: .venv/bin/python plotting/lowest_cusp.py --out plots/
"""
import argparse, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _style as st
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import binom_core as core

EVEN_C, ODD_C = "#1f6fb4", "#d1495b"

def load(path):
    d = np.genfromtxt(path, delimiter=',', names=True)
    n = d['n'].astype(int); i = d['argmin_i'].astype(int); j = d['argmin_j'].astype(int)
    dE = d['min_E_minus_Ehalf']
    p = np.empty(len(n))
    for k in range(len(n)):
        lnC = core.lnC_arr(n[k])
        p[k] = 1.0/(1.0 + np.exp(-(lnC[i[k]] - lnC[j[k]])/(j[k] - i[k])))
    return n, i, j, dE, p, (i + j - n), (j - i)

def frame(a, nrows=1):
    plt.rcParams.update({"font.size": 30, "axes.titlesize": 40, "axes.labelsize": 34,
                         "xtick.labelsize": 26, "ytick.labelsize": 26, "legend.fontsize": 28})
    fig, axes = plt.subplots(nrows, 1, figsize=st.figsize(a), dpi=st.DPI,
                             sharex=True, gridspec_kw=dict(hspace=0.10))
    return fig, (axes if nrows > 1 else [axes])

def setx(ax, a):
    if a.logx: ax.set_xscale('log')

def scatter_parity(ax, n, y, ms, **kw):
    ev = n % 2 == 0
    ax.scatter(n[ev], y[ev], s=ms, c=EVEN_C, alpha=0.75, linewidths=0, rasterized=True, **kw)
    ax.scatter(n[~ev], y[~ev], s=ms, c=ODD_C, alpha=0.75, linewidths=0, rasterized=True, **kw)

def legend(ax, extra=()):
    h = [Line2D([], [], marker='o', ls='', ms=22, color=EVEN_C, label='n even  (width odd)'),
         Line2D([], [], marker='o', ls='', ms=22, color=ODD_C,  label='n odd   (width even)')]
    ax.legend(handles=list(h) + list(extra), loc='best', framealpha=0.92, markerscale=1.0)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", default="public/per_n_summary.csv")
    ap.add_argument("--out", default="plots")
    ap.add_argument("--logx", action="store_true",
                    help="log x axis; default is LINEAR so all ~5000 n are spread evenly")
    st.add_size_args(ap)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    n, i, j, dE, p, band, w = load(a.summary)
    ms = st.marker_size(len(n), width_px=a.width)
    ev = n % 2 == 0
    cE = dE*n**1.5; cP = (p - 0.5)*2*(n + 1)
    m = n >= 1000
    # Report the bands actually present rather than trusting the observation that they are all 1.
    # The parity legend below ("n even -> width odd") is only valid in band 1: w = n + b - 2i, so
    # width parity follows n+b, not n.  If other bands ever appear, that legend needs rethinking.
    bands = np.unique(band)
    bandtxt = (f"band $i+j-n$ = {bands[0]} for all {len(n)} values of $n$"
               if len(bands) == 1 else
               f"bands present: {', '.join(map(str, bands[:6]))}"
               + ("..." if len(bands) > 6 else "") + " -- THE PARITY LEGEND ASSUMES BAND 1")
    sub = (bandtxt + "; "
           f"$(E-E_{{1/2}})\\,n^{{3/2}} \\to$ {np.median(cE[m&ev]):.4f} (even), "
           f"{np.median(cE[m&~ev]):.4f} (odd)")

    # --- 1. E - E(1/2) ------------------------------------------------------
    fig, (ax, ax2) = frame(a, 2)
    scatter_parity(ax, n, dE, ms)
    setx(ax, a); ax.set_yscale('log')
    ax.set_ylabel(r"$E(p^*) - E(1/2)$")
    ax.set_title("Lowest cusp of each $n$: height above $E(1/2)$\n" + sub)
    ax.grid(alpha=0.25, which='both')
    legend(ax)
    scatter_parity(ax2, n, cE, ms)
    setx(ax2, a); ax2.set_xlabel("$n$")
    ax2.set_ylabel(r"$(E(p^*) - E(1/2))\; n^{3/2}$")
    ax2.grid(alpha=0.25, which='both')
    st.save(fig, os.path.join(a.out, "lowest_cusp_E.png")); plt.close(fig)

    # --- 2. p* --------------------------------------------------------------
    fig, (ax, ax2) = frame(a, 2)
    scatter_parity(ax, n, p - 0.5, ms)
    setx(ax, a); ax.set_yscale('log')
    ax.set_ylabel(r"$p^* - 1/2$")
    ax.set_title("Lowest cusp of each $n$: its tie point $p^*$\n"
                 r"every band-1 tie point sits at $p^*-1/2 \approx b/(2(n+1))$, so the whole band "
                 r"is a cluster just above $1/2$")
    ax.grid(alpha=0.25, which='both'); legend(ax)
    scatter_parity(ax2, n, cP, ms)
    ax2.axhline(1.0, color='0.35', lw=2.5, ls='--', zorder=1)
    setx(ax2, a); ax2.set_xlabel("$n$")
    ax2.set_ylabel(r"$(p^* - 1/2)\cdot 2(n{+}1)$")
    ax2.grid(alpha=0.25, which='both')
    st.save(fig, os.path.join(a.out, "lowest_cusp_p.png")); plt.close(fig)

    # --- 3. the colouring quantity: width ------------------------------------
    fig, (ax,) = frame(a, 1)
    scatter_parity(ax, n, w, ms)
    c = np.median(w[m]/np.sqrt(n[m]))
    xs = np.linspace(3, n.max(), 800)
    ax.plot(xs, c*np.sqrt(xs), color='0.25', lw=3, ls='--', zorder=3)
    setx(ax, a)
    if a.logx: ax.set_yscale('log')
    ax.set_xlabel("$n$"); ax.set_ylabel(r"width $j-i$ of the lowest cusp")
    ax.set_title("Lowest cusp of each $n$: its width\n" +
                 ("band is always 1 ($i+j=n+1$), so the width is the only free parameter"
                  if len(bands) == 1 and bands[0] == 1 else
                  f"bands present: {list(bands[:6])}"))
    ax.grid(alpha=0.25, which='both')
    legend(ax, extra=[Line2D([], [], color='0.25', lw=3, ls='--',
                             label=f"${c:.4f}\\,\\sqrt{{n}}$")])
    st.save(fig, os.path.join(a.out, "lowest_cusp_width.png")); plt.close(fig)
    print(f"wrote 3 figures to {a.out}/ at {a.width}x{a.height} "
          f"(marker {np.sqrt(ms)*st.DPI/72:.2f} px vs {a.width*0.88/len(n):.2f} px spacing)")

if __name__ == "__main__":
    main()
