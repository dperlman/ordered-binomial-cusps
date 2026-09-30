"""plotting/E_over_n.py -- how far E sits above E(1/2), per n, at the first two tie points and cusps
above 1/2 and at two fixed p.

    .venv/bin/python plotting/E_over_n.py [--nmax 100] [--absolute] [--out plots]

Plotted: (E(n,p) - E(n,1/2))/n, log y, for n = 2..nmax, at
  first tie point    smallest p* > 1/2  (exhaustive pair scan, first_tie_vs_cusp._two_lowest;
                     E from binom_core.evaluate; from n=2)
  first cusp         smallest-p* cusp   (cusps/cusps_all.csv; from n=3)
  second tie point   second-smallest p* (from n=3)
  second cusp        second-smallest-p* cusp (from n=6: n=3..5 have one cusp above 1/2)
  p = 0.51, 0.61     binom_core.E_at
--series picks which tie/cusp series appear (tK = K-th tie point, cK = K-th cusp) and --p the
fixed p values, e.g.
  --series c1,c2 --p 0.6,0.51,0.501 --tag _cusps
  --series t1,t2,t5,t10,t20,t50,t100 --p "" --tag _ties      (K-th tie point by rank)
--absolute plots E/n itself instead (linear y) -- there the curves overlap almost completely.
"First" excludes p=1/2 itself, where every mirror pair ties at once.  The first tie point is the
first cusp only at n = 3, 4, 5, 6, 7, 9 (first_tie_vs_cusp.py, n<=5000), so those points coincide;
tie points are drawn as "+", cusps as hollow diamonds, fixed p as dots, so coinciding points
stay visible.  All values are positive for
n<=100 (E > E(1/2) everywhere checked).
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
from first_tie_vs_cusp import lowest_k

def cusps_by_n(path, nmax):
    """{n: (E, band) arrays of that n's cusps, in increasing p*}."""
    import pyarrow.csv as pc
    t = pc.read_csv(path, convert_options=pc.ConvertOptions(
        include_columns=["n", "i", "j", "pstar", "E"]))
    n, i, j, p, E = (t[c].to_numpy() for c in ("n", "i", "j", "pstar", "E"))
    m = n <= nmax; n, i, j, p, E = n[m], i[m], j[m], p[m], E[m]
    o = np.lexsort((p, n)); n, E, band = n[o], E[o], (i + j - n)[o]
    ns, first = np.unique(n, return_index=True)
    return {int(a): (e, b) for a, e, b in zip(ns, np.split(E, first[1:]), np.split(band, first[1:]))}

def parse_series(text):
    keys = [k for k in text.split(",") if k]
    if any(k[0] not in "tc" or not k[1:].isdigit() or int(k[1:]) < 1 for k in keys):
        sys.exit("--series entries must look like t1, t20, c2 ...")
    return keys

def series_values(ns, keys, fixed, csv, absolute=False):
    """{key: (n, value, band)} for tie/cusp keys (tK, cK) and fixed p (key f"p{v}", band -1).
    value is (E - E(1/2))/n, or E/n with absolute=True.  Also returns the smallest ln-rho gap
    between consecutive tie ranks 1..K+1 (inf when no tie series was asked for)."""
    K = max([int(k[1:]) for k in keys if k[0] == "t"], default=0)
    cz = cusps_by_n(csv, int(ns[-1]))
    rows = {k: ([], [], []) for k in keys + [f"p{v}" for v in fixed]}
    min_gap = np.inf
    for n in ns:
        n = int(n); base = 0.0 if absolute else core.E_half(n)
        Et, Bt = [], []
        if K:                                          # O(n^2) scan: only when asked for
            lnC = core.lnC_arr(n)
            ti, tj, tv = lowest_k(n, lnC, K)
            if len(tv) > 1: min_gap = min(min_gap, float(np.diff(tv).min()))
            Et = [core.evaluate(n, int(i), int(j), lnC)[1] for i, j in zip(ti[:K], tj[:K])]
            Bt = [int(i + j - n) for i, j in zip(ti[:K], tj[:K])]
        Ec, Bc = cz.get(n, ([], []))
        have = {}
        for k in keys:
            q = int(k[1:]) - 1
            E, B = (Et, Bt) if k[0] == "t" else (Ec, Bc)
            if q < len(E): have[k] = (E[q], B[q])
        have.update({f"p{v}": (core.E_at(n, v), -1) for v in fixed})
        for key, (E, B) in have.items():
            rows[key][0].append(n); rows[key][1].append((E - base)/n); rows[key][2].append(B)
    return {k: tuple(np.array(c) for c in v) for k, v in rows.items()}, min_gap

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="cusps/cusps_all.csv")
    ap.add_argument("--nmin", type=int, default=2)
    ap.add_argument("--nmax", type=int, default=100)
    ap.add_argument("--absolute", action="store_true", help="plot E/n, not (E - E(1/2))/n")
    ap.add_argument("--out", default="plots")
    ap.add_argument("--height", type=int, default=st.DEFAULT_H)
    ap.add_argument("--series", default="t1,c1,t2,c2",
                    help="tie/cusp series: tK = K-th tie point above 1/2, cK = K-th cusp")
    ap.add_argument("--p", default="0.51,0.61", help="fixed p values, drawn as dots")
    ap.add_argument("--tag", default="", help="suffix for the output file name")
    a = ap.parse_args()
    keys = parse_series(a.series)
    K = max([int(k[1:]) for k in keys if k[0] == "t"], default=0)
    fixed = [float(v) for v in a.p.split(",") if v]

    ns = np.arange(a.nmin, a.nmax + 1)
    R, min_gap = series_values(ns, keys, fixed, a.csv, a.absolute)

    # (key, label, colour, size px, shape, stroke px): tie points "+", cusps hollow diamonds,
    # fixed p round dots.  Drawn in this order.
    STYLE = {"t1": ("first tie point", "#2a78d6", 36, "plus", 5),
             "c1": ("first cusp", "#eb6834", 34, "diamond", 4),
             "t2": ("second tie point", "#1baf7a", 26, "plus", 4),
             "c2": ("second cusp", "#4a3aa7", 24, "diamond", 3.5)}
    DOT = {0.51: "#222222", 0.61: "#8a8a8a", 0.6: "#8a8a8a", 0.501: "#1baf7a"}
    spare = iter(["#2a78d6", "#eb6834", "#e87ba4", "#008300"])
    if set(keys) <= set(STYLE):
        series = [(k, *STYLE[k]) for k in keys]
    else:                                              # many ranks: one ordinal ramp per kind
        import matplotlib as mpl
        series = []
        for kind, shape, lw, word in (("t", "plus", 5, "tie point"), ("c", "diamond", 4, "cusp")):
            ks = [k for k in keys if k[0] == kind]
            cols = mpl.colormaps["viridis"](np.linspace(0.0, 0.88, max(len(ks), 2)))
            series += [(k, f"{word} #{k[1:]}", mpl.colors.to_hex(cols[q]), 36, shape, lw)
                       for q, k in enumerate(ks)]
    series += [(f"p{v}", f"p = {v:g}", DOT.get(v) or next(spare), 12, "disc", 0) for v in fixed]

    g = st.NGrid(ns[0], ns[-1], height=a.height)
    # marks and strokes were sized for 60 px per n; shrink with the column, with floors
    sc = min(1.0, max(0.3, g.k/60))
    series = [(k, lab, col, max(4, round(d*sc)), shape, max(1.3, lw*sc))
              for k, lab, col, d, shape, lw in series]
    ax, b = g.ax, g.back[0]
    allv = np.concatenate([R[k][1] for k, *_ in series])
    if a.absolute:
        pad = 0.03*(allv.max() - allv.min()); ax.set_ylim(allv.min() - pad, allv.max() + pad)
    else:
        if (allv <= 0).any(): sys.exit("a value is <= 0: E fell below E(1/2) -- check before plotting")
        ax.set_yscale("log"); ax.set_ylim(allv.min()/1.6, allv.max()*1.6)
    for key, _, col, d, shape, lw in series:
        x, y, _ = R[key]
        b.plot(x, y, color=col, lw=2 if g.k >= 20 else 1.1, alpha=0.55, zorder=2)   # under marks
        g.marks(x, y, col, d, shape, lw=lw)

    what = "$E/n$" if a.absolute else "$(E(n,p) - E(n,\\frac{1}{2}))\\,/\\,n$"
    parts = {"t1": "first tie point", "t2": "second tie point", "c1": "first cusp",
             "c2": "second cusp"}
    if set(keys) <= set(parts):
        head = ", ".join(parts[k] for k in keys)
    else:
        grp = [(w, [k[1:] for k in keys if k[0] == c]) for c, w in (("t", "tie points"), ("c", "cusps"))]
        head = "; ".join(f"{w} #{', '.join(v)}" for w, v in grp if v)
    tail = "$p = " + ",\\ ".join(f"{v:g}" for v in fixed) + "$" if fixed else ""
    b.set_title(f"{what}  at the {head}" + (" (above $\\frac{1}{2}$)" if keys else "")
                + (f" and at {tail}" if tail else ""), fontsize=44, pad=30)
    b.set_xlabel("$n$", fontsize=38, labelpad=18)
    b.set_ylabel(what, fontsize=38, labelpad=22)
    b.tick_params(labelsize=28, length=12, width=2)
    b.tick_params(which="minor", length=6, width=1.2)
    from matplotlib.ticker import MultipleLocator
    span = ns[-1] - ns[0]
    major, minor = ((10, 1) if span <= 150 else (50, 10) if span <= 600 else
                    (100, 20) if span <= 1500 else (500, 100))
    b.xaxis.set_major_locator(MultipleLocator(major)); b.xaxis.set_minor_locator(MultipleLocator(minor))
    b.grid(True, which="major", alpha=0.35, lw=1.2)
    if not a.absolute: b.grid(True, which="minor", axis="y", alpha=0.15, lw=0.8)
    b.set_axisbelow(True)
    mk = {"plus": "+", "diamond": "D", "disc": "o"}
    handles = [Line2D([], [], marker=mk[shape], color=col, lw=2, alpha=0.9,
                      ms=(d if shape != "diamond" else d/np.sqrt(2))*72/st.DPI,
                      mew=lw*72/st.DPI if lw else 1, mfc="none" if shape == "diamond" else col,
                      label=f"{lab}  (n = {R[k][0][0]}..{R[k][0][-1]})")
               for k, lab, col, d, shape, lw in series]
    ax.legend(handles=handles, fontsize=30, loc="lower left" if not a.absolute else "lower right",
              framealpha=0.93)

    os.makedirs(a.out, exist_ok=True)
    name = f"E_{'over_n' if a.absolute else 'minus_half_over_n'}_n{a.nmax:05d}{a.tag}.png"
    path = g.save(os.path.join(a.out, name))
    print(f"{path}  {g.W}x{g.H} px, {g.k} px per n")
    if K: print(f"  smallest ln-rho gap between consecutive tie ranks 1..{K+1}: {min_gap:.3e}")
    for k, lab, *_ in series:
        x, y, _ = R[k]
        print(f"  {lab:<17} n={x[0]}..{x[-1]}  {y[0]:.4e} .. {y[-1]:.4e}")

if __name__ == "__main__":
    main()
