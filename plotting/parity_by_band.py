"""plotting/parity_by_band.py -- one tie rank's E and its even/odd offset, every point coloured by
the BAND (i+j-n) of the pair that holds that rank at that n.

    .venv/bin/python plotting/parity_by_band.py [--rank 100] [--nmin 20] [--nmax 400]

Top:    y(n) = (E(n,p*) - E(n,1/2))/n at the rank-K tie point above 1/2 (log y).
Bottom: the band itself as a staircase, so higher and lower bands read directly.
Middle: alpha(n) = ((-1)^n/2)[ln y(n) - (ln y(n-1) + ln y(n+1))/2], with NO points left out --
        at a band switch it measures the jump, which is part of what is being shown.
A band-b pair has width n+b-2i, so its width parity is the parity of n+b: within one band, even and
odd n always get pairs of opposite width parity (the period-2 lock).  Which parity comes out higher
depends on the band (RESEARCH_LOG 2026-09-29), and which band holds rank K changes as the band
families slide past each other -- this plot shows both at once.
Bands 1-4 get a validated categorical palette; 5 and above are grouped in grey.  Each n's column is
also tinted in its band's colour, so the band regions read as vertical stripes behind both panels.
"""
import argparse, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _style as st
from E_over_n import series_values
from parity_alternation import parity_offset

BAND_COL = {1: "#2a78d6", 2: "#eb6834", 3: "#1baf7a", 4: "#4a3aa7"}
OTHER = "#8a8a8a"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="cusps/cusps_all.csv")
    ap.add_argument("--rank", type=int, default=100)
    ap.add_argument("--nmin", type=int, default=20)
    ap.add_argument("--nmax", type=int, default=400)
    ap.add_argument("--out", default="plots")
    ap.add_argument("--height", type=int, default=st.DEFAULT_H)
    a = ap.parse_args()

    key = f"t{a.rank}"
    R, _ = series_values(np.arange(2, a.nmax + 2), [key], [], a.csv)   # +1 so alpha(nmax) exists
    x, y, band = R[key]
    al, sw = parity_offset(x, y, band, 2, keep=True)
    m = (x >= a.nmin) & (x <= a.nmax)
    x, y, band, al = x[m], y[m], band[m], al[m]
    col = np.array([BAND_COL.get(int(b), OTHER) for b in band])

    g = st.NGrid(int(x[0]), int(x[-1]), height=a.height, nrows=3, gap=170)
    top, bot, bnd = g.axes
    bnd.set_ylim(0.4, band.max() + 0.6)
    top.set_yscale("log"); top.set_ylim(y.min()/1.5, y.max()*1.5)
    ok = np.isfinite(al)
    lo, hi = al[ok].min(), al[ok].max(); pad = 0.06*(hi - lo)
    bot.set_ylim(lo - pad, hi + pad)
    sc = min(1.0, max(0.3, g.k/60)); d = max(7, round(40*sc)); lw = max(1.6, 5*sc)
    for panel in (0, 1, 2):                               # background: each n's column in its band tint
        for c in np.unique(col):
            g.column_fill(x[col == c], c, panel=panel, alpha=0.16)
    for panel, v in ((0, y), (1, al)):
        b = g.back[panel]
        if panel == 1: b.axhline(0.0, color="0.3", lw=2, zorder=1.5)
        b.plot(x, v, color="0.75", lw=1.3, zorder=2)                     # under the marks
        for c in np.unique(col):
            q = (col == c) & np.isfinite(v)
            g.marks(x[q], v[q], c, d, "plus", lw=lw, panel=panel)
    # the band itself as a staircase: steps at n +- 1/2, i.e. exactly on the column edges
    g.back[2].plot(x, band, drawstyle="steps-mid", color="0.35", lw=2.2, zorder=2)
    for c in np.unique(col):
        q = col == c
        g.marks(x[q], band[q].astype(float), c, max(6, round(d*0.6)), "disc", panel=2)

    b0, b1, b2 = g.back
    b0.set_title(f"Tie point #{a.rank} above $\\frac{{1}}{{2}}$, coloured by the band $i+j-n$ of its pair",
                 fontsize=42, pad=28)
    b0.set_ylabel("$(E - E(\\frac{1}{2}))\\,/\\,n$", fontsize=34, labelpad=18)
    b1.set_title("even/odd offset  $\\alpha(n) = \\frac{(-1)^n}{2}\\,[\\,\\ln y_n - \\frac{1}{2}"
                 "(\\ln y_{n-1} + \\ln y_{n+1})\\,]$,  no points left out", fontsize=36, pad=20)
    b1.set_ylabel("$\\alpha$   (> 0: even $n$ higher)", fontsize=32, labelpad=18)
    b2.set_title("the band $i+j-n$ of the pair holding that rank", fontsize=36, pad=20)
    b2.set_ylabel("band", fontsize=32, labelpad=18)
    b2.set_xlabel("$n$", fontsize=38, labelpad=18)
    from matplotlib.ticker import MaxNLocator
    b2.yaxis.set_major_locator(MaxNLocator(integer=True))
    from matplotlib.ticker import MultipleLocator
    span = x[-1] - x[0]
    major, minor = ((10, 1) if span <= 150 else (50, 10) if span <= 600 else
                    (100, 20) if span <= 1500 else (500, 100))
    for b in g.back:
        b.tick_params(labelsize=28, length=12, width=2)
        b.tick_params(which="minor", length=6, width=1.2)
        b.xaxis.set_major_locator(MultipleLocator(major)); b.xaxis.set_minor_locator(MultipleLocator(minor))
        b.grid(True, alpha=0.35, lw=1.2); b.set_axisbelow(True)
    counts = {c: int((col == c).sum()) for c in set(col)}
    labels = [(BAND_COL[k], f"band {k}") for k in sorted(BAND_COL)] + [(OTHER, "band 5 and above")]
    from matplotlib.patches import Patch
    handles = [(Patch(fc=c, alpha=0.3, ec="none"),
                Line2D([], [], marker="+", ls="", color=c, ms=18, mew=3.5)) for c, lab in labels
               if counts.get(c)]
    names = [f"{lab}  ({counts.get(c, 0)} n)" for c, lab in labels if counts.get(c)]
    top.legend(handles, names, fontsize=28, loc="upper right", framealpha=0.93, ncol=len(handles),
               handlelength=2.2)

    os.makedirs(a.out, exist_ok=True)
    path = g.save(os.path.join(a.out, f"parity_by_band_rank{a.rank:03d}_n{a.nmin}-{a.nmax}.png"))
    print(f"{path}  {g.W}x{g.H} px, {g.k} px per n")
    print(f"  bands: {dict(zip(*np.unique(band, return_counts=True)))};  band switches in range: "
          f"{sum(1 for n in sw if a.nmin <= n <= a.nmax)}")

if __name__ == "__main__":
    main()
