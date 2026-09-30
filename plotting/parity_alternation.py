"""plotting/parity_alternation.py -- the even/odd alternation in (E - E(1/2))/n, demodulated.

    .venv/bin/python plotting/parity_alternation.py [--nmax 1000] [--series t1,t2,t5,...] [--p ""]

Takes the same series as E_over_n.py (series_values) and, for each, the parity offset

    alpha(n) = ((-1)^n / 2) * [ ln y(n) - (ln y(n-1) + ln y(n+1))/2 ],    y = (E - E(1/2))/n.

Why this form.  The second difference y(n+1) - 2y(n) + y(n-1) has gain 4 sin^2(w/2) on a
component of angular frequency w: maximal (4) at period 2, near zero for smooth variation.  So it
isolates the even/odd alternation -- but it also passes anything abrupt.  Working in ln y makes the
result a RATIO (independent of the seven decades y spans), and multiplying by (-1)^n turns the
alternating +-4a into a smooth curve.  For y = s(n) exp((-1)^n alpha) it returns alpha up to
O(s''/s).  alpha > 0: even n sit above odd n; alpha < 0: odd n above even.

--order 1 uses the FIRST difference instead: alpha(n) = ((-1)^n/2)[ln y(n) - ln y(n-1)].  For
y = s exp((-1)^n alpha) that is alpha + ((-1)^n/2) Delta ln s: the smooth trend is NOT cancelled but
leaks through as a residual zigzag of size ~|d ln y/dn|/2.  (A ratio panel, e^(2 alpha), was
dropped 2026-09-28: it is a monotone rescaling of alpha and showed nothing new.)

Gaps.  alpha(n) is left out wherever the series' pair changes BAND between n-1 and n+1 (a family
switch, e.g. tie rank 50 at n=104): there the second difference measures the jump, not parity.
Those n are printed.
"""
import argparse, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _style as st
from E_over_n import series_values, parse_series

def parity_offset(x, y, band, order=2, keep=False):
    """alpha on the n that have the needed neighbours; NaN where the band changes.
    order=2: ((-1)^n/2)[ln y(n) - (ln y(n-1) + ln y(n+1))/2]   (n-1, n, n+1)
    order=1: ((-1)^n/2)[ln y(n) - ln y(n-1)]                  (n-1, n)"""
    pos = {int(n): q for q, n in enumerate(x)}
    al = np.full(len(x), np.nan); sw = []
    ly = np.log(y)
    for q, n in enumerate(x):
        a = pos.get(int(n) - 1)
        b = pos.get(int(n) + 1) if order == 2 else q
        if a is None or b is None: continue
        if not (band[a] == band[q] == band[b]):
            sw.append(int(n))
            if not keep: continue
        d = ly[q] - (0.5*(ly[a] + ly[b]) if order == 2 else ly[a])
        al[q] = ((-1)**int(n))/2*d
    return al, sw

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="cusps/cusps_all.csv")
    ap.add_argument("--nmin", type=int, default=2)
    ap.add_argument("--nmax", type=int, default=1000)
    ap.add_argument("--series", default="t1,t2,t5,t10,t20,t50,t100")
    ap.add_argument("--p", default="", help="fixed p values (optional)")
    ap.add_argument("--from-n", type=int, default=20, help="first n shown (small n swamp the scale)")
    ap.add_argument("--out", default="plots")
    ap.add_argument("--height", type=int, default=st.DEFAULT_H)
    ap.add_argument("--tag", default="")
    ap.add_argument("--keep-switches", action="store_true",
                    help="also plot the n where the pair changes band (normally left out)")
    ap.add_argument("--order", type=int, choices=(1, 2), default=2,
                    help="2: second difference (smooth trend cancels); 1: first difference")
    a = ap.parse_args()
    keys = parse_series(a.series)
    fixed = [float(v) for v in a.p.split(",") if v]

    ns = np.arange(a.nmin, a.nmax + 1)
    R, _ = series_values(ns, keys, fixed, a.csv)
    import matplotlib as mpl
    ties = [k for k in keys if k[0] == "t"]; cus = [k for k in keys if k[0] == "c"]
    style = {}
    for group, shape, word in ((ties, "plus", "tie point"), (cus, "diamond", "cusp")):
        cols = mpl.colormaps["viridis"](np.linspace(0.0, 0.88, max(len(group), 2)))
        for q, k in enumerate(group):
            style[k] = (f"{word} #{k[1:]}", mpl.colors.to_hex(cols[q]), shape)
    for v, col in zip(fixed, ["#222222", "#8a8a8a", "#e87ba4", "#008300"]):
        style[f"p{v}"] = (f"p = {v:g}", col, "disc")

    A = {}
    for k in keys + [f"p{v}" for v in fixed]:
        x, y, band = R[k]
        al, sw = parity_offset(x, y, band, a.order, a.keep_switches)
        m = x >= a.from_n
        A[k] = (x[m], al[m], [n for n in sw if n >= a.from_n])

    g = st.NGrid(a.from_n, ns[-1], height=a.height)
    top, b0 = g.ax, g.back[0]
    allv = np.concatenate([v[1][np.isfinite(v[1])] for v in A.values()])
    lo, hi = np.percentile(allv, [0.2, 99.8]); pad = 0.08*(hi - lo)
    top.set_ylim(min(lo, 0) - pad, max(hi, 0) + pad)
    sc = min(1.0, max(0.3, g.k/60)); d = max(5, round(34*sc)); lw = max(1.4, 4.5*sc)
    b0.axhline(0.0, color="0.3", lw=2, zorder=1.5)
    for k, (x, al, _) in A.items():
        lab, col, shape = style[k]
        b0.plot(x, al, color=col, lw=1.4, alpha=0.6, zorder=2)
        ok = np.isfinite(al)
        g.marks(x[ok], al[ok], col, d if shape != "disc" else max(4, round(12*sc)), shape, lw=lw)

    formula = ("[\\,\\ln y_n - \\frac{1}{2}(\\ln y_{n-1} + \\ln y_{n+1})\\,]" if a.order == 2
               else "[\\,\\ln y_n - \\ln y_{n-1}\\,]")
    b0.set_title("Even/odd alternation of $(E - E(\\frac{1}{2}))/n$"
                 + (" (first difference)" if a.order == 1 else "")
                 + f":  $\\alpha(n) = \\frac{{(-1)^n}}{{2}}\\,{formula}$", fontsize=42, pad=28)
    b0.set_ylabel("$\\alpha$   (> 0: even $n$ higher)", fontsize=34, labelpad=18)
    b0.set_xlabel("$n$", fontsize=38, labelpad=18)
    from matplotlib.ticker import MultipleLocator
    span = ns[-1] - a.from_n
    major, minor = ((10, 1) if span <= 150 else (50, 10) if span <= 600 else
                    (100, 20) if span <= 1500 else (500, 100))
    for b in g.back:
        b.tick_params(labelsize=28, length=12, width=2)
        b.tick_params(which="minor", length=6, width=1.2)
        b.xaxis.set_major_locator(MultipleLocator(major)); b.xaxis.set_minor_locator(MultipleLocator(minor))
        b.grid(True, alpha=0.35, lw=1.2); b.set_axisbelow(True)
    mk = {"plus": "+", "diamond": "D", "disc": "o"}
    handles = [Line2D([], [], marker=mk[style[k][2]], color=style[k][1], lw=1.6, ms=14, mew=2.5,
                      mfc="none" if style[k][2] == "diamond" else style[k][1], label=style[k][0])
               for k in A]
    top.legend(handles=handles, fontsize=26, loc="upper right", framealpha=0.93, ncol=len(A))
    nsw = sum(len(v[2]) for v in A.values())
    where = f"between n-1 and n{'+1' if a.order == 2 else ''}"
    top.text(0.005, 0.02, (f"NO gaps: includes {nsw} points where the pair changes band {where}"
                           if a.keep_switches else
                           f"gaps: the pair changes band {where} ({nsw} points left out)"),
             transform=top.transAxes, fontsize=24, color="0.35", va="bottom")

    os.makedirs(a.out, exist_ok=True)
    path = g.save(os.path.join(a.out, f"parity_alternation_n{ns[-1]:05d}"
                               f"{'_d1' if a.order == 1 else ''}"
                               f"{'_nogaps' if a.keep_switches else ''}{a.tag}.png"))
    print(f"{path}  {g.W}x{g.H} px, {g.k} px per n")
    for k, (x, al, sw) in A.items():
        ok = np.isfinite(al)
        tail = al[ok][x[ok] >= ns[-1] - 100]
        print(f"  {style[k][0]:<16} alpha over the last 100 n: even {np.nanmean(tail[x[ok][x[ok] >= ns[-1]-100] % 2 == 0]):+.4f}"
              f" / odd {np.nanmean(tail[x[ok][x[ok] >= ns[-1]-100] % 2 == 1]):+.4f}   band switches at {sw[:8]}{' ...' if len(sw) > 8 else ''}")

if __name__ == "__main__":
    main()
