"""plotting/max_pstar.py -- the largest cusp p* of each n, versus n (linear, one pixel column per n).

    .venv/bin/python plotting/max_pstar.py [--csv cusps/cusps_all.csv] [--data data] [--out plots]

Two panels on the same n columns: the full range (top) and a zoom on the band the maximum settles
into (bottom).  Each n's mark is coloured by the width j-i of the pair that produces its maximum:
nearly always a narrow tie (width 1-4), occasionally a much wider one (width >= 57, first at
n=366).  The limit is fitted as p*_max ~ L - c n^-a over n >= --fit-from, and n=6000/7000/8000 (from
the Parquet cusp dumps, outside the complete table) are printed as an out-of-sample check.

Bears on CLAUDE.md's open item "no cusps above p ~ 0.66" (data, not proof): the maximum is highest
at small n (0.65693 at n=15, pair (9,11)); for n >= 500 it sits in a band that closes on ~0.6522
from both sides.
"""
import argparse, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _style as st
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# width classes: 1, 2, 3-4, wide.  Validated as a set (all pairs, light surface).
CLASSES = [("width 1", "#2a78d6"), ("width 2", "#eb6834"), ("width 3-4", "#1baf7a"),
           ("wide (width >= 5)", "#4a3aa7")]

def load(path):
    """Per n: (n, i, j, p*) of the cusp with the largest p*."""
    import pyarrow.csv as pc
    t = pc.read_csv(path, convert_options=pc.ConvertOptions(include_columns=["n", "i", "j", "pstar"]))
    n, i, j, p = (t[c].to_numpy() for c in ("n", "i", "j", "pstar"))
    o = np.lexsort((-p, n)); n, i, j, p = n[o], i[o], j[o], p[o]
    ns, first = np.unique(n, return_index=True)
    return ns, i[first], j[first], p[first]

def width_class(w):
    return np.select([w == 1, w == 2, w <= 4], [0, 1, 2], 3)

def fit_limit(ns, P, lo):
    """Least squares P ~ L - c n^-a over n >= lo, a on a grid.  Returns L, c, a, rms."""
    m = ns >= lo; best = None
    for a in np.linspace(0.2, 2.0, 181):
        X = np.column_stack([np.ones(m.sum()), -ns[m]**-a])
        co = np.linalg.lstsq(X, P[m], rcond=None)[0]
        rs = ((X @ co - P[m])**2).mean()
        if best is None or rs < best[0]: best = (rs, co[0], co[1], a)
    rs, L, c, a = best
    return L, c, a, np.sqrt(rs)

def parquet_max(data, ns):
    import cusps_data as cd
    out = []
    for n in ns:
        try: d = cd.load_cusps(n, data=data)
        except FileNotFoundError: continue
        p = d["pstar"][~d["is_axis"]]; out.append((n, p.max()))
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="cusps/cusps_all.csv")
    ap.add_argument("--data", default="data")
    ap.add_argument("--out", default="plots")
    ap.add_argument("--height", type=int, default=st.DEFAULT_H)
    ap.add_argument("--fit-from", type=int, default=500)
    a = ap.parse_args()

    ns, I, J, P = load(a.csv)
    if not np.array_equal(ns, np.arange(ns[0], ns[-1] + 1)):
        sys.exit("the table has gaps in n; this plot is for a complete range")
    cls = width_class(J - I)
    L, c, ex, rms = fit_limit(ns, P, a.fit_from)
    ext = parquet_max(a.data, [n for n in (6000, 7000, 8000) if n > ns[-1]])
    kmax = np.argmax(P)

    g = st.NGrid(ns[0], ns[-1], height=a.height, nrows=2, gap=190)
    top, bot = g.axes
    top.set_ylim(0.58, 0.66)
    zlo, zhi = 0.6505, 0.6530
    bot.set_ylim(zlo, zhi)
    for panel, h in ((0, 5), (1, 5)):
        for q, (_, col) in enumerate(CLASSES):
            m = cls == q
            g.points(ns[m], P[m], col, panel=panel, h=h)

    xf = np.arange(a.fit_from, ns[-1] + 1, dtype=float)
    for ax in (top, bot):
        ax.axhline(L, color="0.25", lw=2.5, ls="--")
        ax.plot(xf, L - c*xf**-ex, color="black", lw=3, ls=":")
    top.annotate(f"largest cusp p* at any n:  {P[kmax]:.5f}  (n = {ns[kmax]}, pair "
                 f"({I[kmax]},{J[kmax]}))", xy=(ns[kmax], P[kmax]),
                 xytext=(ns[kmax] + 0.06*len(ns), 0.638), fontsize=28,
                 arrowprops=dict(arrowstyle="->", lw=2.5, color="0.3"))
    nwide = int((cls == 3).sum()); fw = ns[cls == 3][0] if nwide else None
    spread = [(lo, P[(ns >= lo) & (ns < lo + 500)]) for lo in (1000, ns[-1] - 499)]
    spread = [(lo, v.min(), v.max()) for lo, v in spread]
    chk = "   ".join(f"n={n}: {p:.5f} (fit {L - c*n**-ex:.5f})" for n, p in ext)
    top.text(0.985, 0.62,
             f"Fit over n $\\geq$ {a.fit_from}:  $p^*_{{max}} \\approx {L:.5f} - {c:.3f}\\,n^{{-{ex:.2f}}}$"
             f"   (rms {rms:.1e})\n"
             f"Out of sample (Parquet): {chk}\n"
             f"The band closes on L from both sides:  n = {spread[0][0]}-{spread[0][0]+499}: "
             f"{spread[0][1]:.5f}-{spread[0][2]:.5f},   n = {spread[1][0]}-{spread[1][0]+499}: "
             f"{spread[1][1]:.5f}-{spread[1][2]:.5f}.\n"
             f"The maximum comes from a narrow tie at {len(ns) - nwide} of {len(ns)} n; "
             f"from a wide one (width {int((J-I)[cls == 3].min())}-{int((J-I)[cls == 3].max())}) "
             f"at {nwide}, first at n = {fw}.",
             transform=top.transAxes, ha="right", va="top", fontsize=26, color="0.2",
             bbox=dict(boxstyle="round,pad=0.6", fc="white", ec="0.7", alpha=0.93))

    b0, b1 = g.back
    b0.set_title("The largest cusp $p^*$ of each $n$", fontsize=44, pad=30)
    b1.set_title(f"zoom:  {zlo} $\\leq p^* \\leq$ {zhi}", fontsize=34, pad=18)
    b1.set_xlabel("$n$", fontsize=38, labelpad=18)
    for b in g.back:
        b.set_ylabel("max $p^*$ over cusps", fontsize=34, labelpad=22)
        b.tick_params(labelsize=28, length=12, width=2)
        b.grid(True, alpha=0.35, lw=1.2); b.set_axisbelow(True)

    mark = lambda col, l: Line2D([], [], marker="s", ls="", ms=13, mfc=col, mec=col, label=l)
    counts = np.bincount(cls, minlength=4)
    handles = [mark(col, f"{lab}  ({counts[q]:,} n)") for q, (lab, col) in enumerate(CLASSES)]
    handles += [Line2D([], [], lw=3, ls=":", color="black", label="fit  $L - c\\,n^{-a}$"),
                Line2D([], [], lw=2.5, ls="--", color="0.25", label=f"fitted limit L = {L:.5f}")]
    top.legend(handles=handles, fontsize=27, loc="lower right", framealpha=0.93, ncol=2)

    os.makedirs(a.out, exist_ok=True)
    path = g.save(os.path.join(a.out, "max_pstar.png"))
    print(f"{path}  {g.W}x{g.H} px, {g.k} px per n  ({os.path.getsize(path)/1e6:.1f} MB)")
    print(f"  n = {ns[0]}..{ns[-1]};  largest {P[kmax]:.6f} at n={ns[kmax]} ({I[kmax]},{J[kmax]})")
    print(f"  fit n>={a.fit_from}: L={L:.6f} c={c:.4g} a={ex:.2f} rms={rms:.2e}")
    print(f"  width classes {dict(zip([l for l, _ in CLASSES], counts.tolist()))}")
    print(f"  outside zoom (n>={a.fit_from}): {int(((P < zlo) | (P > zhi))[ns >= a.fit_from].sum())}")
    for n, p in ext: print(f"  Parquet n={n}: {p:.6f}  fit {L - c*n**-ex:.6f}")
    for lo, a_, b_ in spread: print(f"  n={lo}..{lo+499}: {a_:.6f} .. {b_:.6f}  (spread {b_-a_:.2e})")

if __name__ == "__main__":
    main()
