"""plotting/E_over_n.py -- E/n against n for a few fixed choices of p, on one shared scale.

    .venv/bin/python plotting/E_over_n.py [--nmax 100] [--out plots]

Series, for n = 2..nmax, counting tie points and cusps upward from p = 1/2:
  first tie = first cusp   p = 1/2 itself: every mirror pair (i,n-i) ties there, and it is always a
                           cusp (binom_core.axis_point).  E = E_half.
  second tie point         the tie point with the smallest p* > 1/2 (binom_core.screen)
  second cusp              the cusp with the smallest p* > 1/2 (cusps/cusps_all.csv; none at n=2)
  p = 0.51, 0.61           binom_core.E_at
E at the tie points is the kernel's (normalised masses); E_at uses the same convention and agrees
with it to ~1e-14.  The second tie point IS the second cusp only at n = 3, 4, 5, 6, 7, 9 (n<=100);
the series are drawn largest-first so coinciding points stay visible.
"""
import argparse, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _style as st
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import binom_core as core

def lowest_cusps(path, nmax):
    """{n: (i, j, p*, E)} for the smallest-p* cusp of each n <= nmax."""
    import pyarrow.csv as pc
    t = pc.read_csv(path, convert_options=pc.ConvertOptions(
        include_columns=["n", "i", "j", "pstar", "E"]))
    n, i, j, p, E = (t[c].to_numpy() for c in ("n", "i", "j", "pstar", "E"))
    m = n <= nmax; n, i, j, p, E = n[m], i[m], j[m], p[m], E[m]
    o = np.lexsort((p, n)); n, i, j, p, E = n[o], i[o], j[o], p[o], E[o]
    ns, first = np.unique(n, return_index=True)
    return {int(a): (int(i[k]), int(j[k]), p[k], E[k]) for a, k in zip(ns, first)}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="cusps/cusps_all.csv")
    ap.add_argument("--nmin", type=int, default=2)
    ap.add_argument("--nmax", type=int, default=100)
    ap.add_argument("--out", default="plots")
    ap.add_argument("--height", type=int, default=st.DEFAULT_H)
    a = ap.parse_args()

    ns = np.arange(a.nmin, a.nmax + 1)
    tie = np.array([(lambda s: s["E"][np.argmin(s["pstar"])])(core.screen(int(n), collect_all=True))
                    for n in ns])/ns
    half = np.array([core.E_half(int(n)) for n in ns])/ns
    fc = lowest_cusps(a.csv, a.nmax)
    cn = np.array([n for n in ns if n in fc])
    cusp = np.array([fc[n][3] for n in cn])/cn
    e51 = np.array([core.E_at(int(n), 0.51) for n in ns])/ns
    e61 = np.array([core.E_at(int(n), 0.61) for n in ns])/ns

    # (label, n, E/n, colour, circle diameter in px) -- drawn in this order, largest first
    series = [("second tie point", ns, tie, "#2a78d6", 30),
              ("second cusp", cn, cusp, "#eb6834", 20),
              ("first tie = first cusp  (p = 1/2)", ns, half, "#222222", 15),
              ("p = 0.51", ns, e51, "#1baf7a", 12),
              ("p = 0.61", ns, e61, "#4a3aa7", 12)]

    g = st.NGrid(ns[0], ns[-1], height=a.height)
    ax, b = g.ax, g.back[0]
    allv = np.concatenate([s[2] for s in series])
    pad = 0.03*(allv.max() - allv.min())
    ax.set_ylim(allv.min() - pad, allv.max() + pad)
    for _, x, y, col, d in series:
        b.plot(x, y, color=col, lw=2, alpha=0.55, zorder=2)      # back axes: under the circles
        g.discs(x, y, col, d)

    b.set_title("$E(n,p)/n$ at the first and second tie points and cusps, $p = 0.51$ and $p = 0.61$",
                fontsize=44, pad=30)
    b.set_xlabel("$n$", fontsize=38, labelpad=18)
    b.set_ylabel("$E/n$", fontsize=38, labelpad=22)
    b.tick_params(labelsize=28, length=12, width=2)
    b.tick_params(which="minor", length=6, width=1.2)
    from matplotlib.ticker import MultipleLocator
    b.xaxis.set_major_locator(MultipleLocator(10)); b.xaxis.set_minor_locator(MultipleLocator(1))
    b.grid(True, alpha=0.35, lw=1.2); b.set_axisbelow(True)
    handles = [Line2D([], [], marker="o", color=col, lw=2, alpha=0.9, ms=d*72/st.DPI, label=lab)
               for lab, _, _, col, d in series]
    ax.legend(handles=handles, fontsize=30, loc="lower right", framealpha=0.93)

    os.makedirs(a.out, exist_ok=True)
    path = g.save(os.path.join(a.out, f"E_over_n_n{a.nmax:05d}.png"))
    print(f"{path}  {g.W}x{g.H} px, {g.k} px per n")
    for lab, x, y, _, _ in series:
        print(f"  {lab:<34} n={x[0]}..{x[-1]}  E/n {y[0]:.5f} .. {y[-1]:.5f}")

if __name__ == "__main__":
    main()
