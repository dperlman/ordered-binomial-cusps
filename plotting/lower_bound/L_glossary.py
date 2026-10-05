"""plotting/lower_bound/L_glossary.py -- a labelled picture of the terms used for L, for explanation only.

    .venv/bin/python plotting/lower_bound/L_glossary.py [--n 12] [--out DIR]

One small n, a window just above p = 1/2 holding two switch points (one HALF, one INT), three
cusps and the certificates that take turns being active.  Every object in the window is labelled
and defined in the glossary beside it.  Left: the whole of [0,1] with the window boxed, and the
glossary.  Right: the window, values relative to E(1/2).  The annotation positions are tuned for
n = 12; other n will draw correctly but the labels may need moving.
"""
import argparse, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _data import OUT, core, switch_points
from L_small_n import tie_points, cusps, family

CE, CL, CS = "#111111", "#d1495b", "#9aa5b1"
CC = {5.75: "#7b6fb0", 6.25: "#2a9d8f", 6.75: "#e9a03b", 7.25: "#3a7dc9", 7.75: "#8c8c8c"}

GLOSSARY = r"""$\bf{Glossary}$   (n fixed;  $f_p(k) = \binom{n}{k} p^k (1-p)^{n-k}$,  $K \sim \mathrm{Bin}(n,p)$)

$\bf{E(n,p)}$   the ordered-binomial expectation $\sum_k w_k\, f_p(k)$, with $w_k$ the
        rank of $f_p(k)$ among the $n+1$ masses (0 = smallest).  Black.

$\bf{tie\ point}\ (i,j)$   the $p^*$ where $f_p(i) = f_p(j)$: two masses swap rank, and
        E gets a convex kink.  Ticks on E.

$\bf{cusp}$   a tie point that is a local minimum of E.  Open circles.

$\bf{certificate}\ W_c$   E scored with a FIXED ranking: $k$ ranked by its distance
        from a centre $c \in \mathbb{Z}/2 + 1/4$ (nearest = rank $n$).  A fixed ranking
        can never beat the true one, so $W_c \leq E$ at every $p$.
        $W_c(p) = n + 1/2 - 2\,E_p|K - c|$.   Thin coloured curves.

$\bf{L}$  (the $\bf{certificate\ envelope}$)   $L(p) = \max_c W_c(p) = n + 1/2 - 2D(p)$,
        $D(p) = \min_c E_p|K - c|$.   $L \leq E$ everywhere, $L(1/2) = E(1/2)$.  Red.

$\bf{active\ certificate}$   the $W_c$ that equals L on a given stretch of $p$.

$\bf{switch\ point}\ p_m$   where the active centre moves up by 1/2; a kink of L, and
        a local minimum.  m-th one at $x \approx m/2$.  Two types, alternating:
        HALF  $P(K \leq a) = 1/2$,   INT  $P(K < a) = P(K > a)$.

$\bf{slack}\ E - L \geq 0$   zero where the distance ranking IS the true ranking.  Grey.

$\bf{margin\ of\ L}$   $L(p) - L(1/2)$.   ($\bigstar$): margin $\geq 0$ for all $p$
        $\Rightarrow$ $E(p) \geq E(1/2)$ for all $p$.

$\bf{x}$ $= n(p - 1/2)$, the scaled distance from 1/2 used in the other L plots."""

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=12)
    ap.add_argument("--lo", type=float, default=0.4965)
    ap.add_argument("--hi", type=float, default=0.5925)
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    n, lo, hi = a.n, a.lo, a.hi

    # ---- data
    m, aa, half, sw, _ = switch_points(n)
    ties = tie_points(n)
    l = core.lnC_arr(n)
    I, J = np.triu_indices(n + 1, 1)
    pairs = {round(float(v), 12): (int(i), int(j)) for i, j, v in zip(I, J, ties)}
    cp, _ = cusps(n)
    p = np.unique(np.concatenate([np.linspace(lo, hi, 4000), ties[(ties >= lo) & (ties <= hi)],
                                  sw[(sw >= lo) & (sw <= hi)], [0.5]]))
    E = np.array([core.E_at(n, v) for v in p])
    c, W = family(n, p)
    L = W.max(axis=0)
    Eh = core.E_half(n)
    yE, yL, yW = E - Eh, L - Eh, W - Eh
    at = lambda arr, v: float(np.interp(v, p, arr))

    pf = np.unique(np.concatenate([np.linspace(0, 1, 3001)[1:-1], ties, sw, 1 - sw]))
    Ef = np.array([core.E_at(n, v) for v in pf])
    Lf = family(n, pf)[1].max(axis=0)

    # ---- figure
    plt.rcParams.update({"font.size": 19})
    fig = plt.figure(figsize=(34, 19), dpi=100)
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.75], height_ratios=[1, 2.1], wspace=0.08, hspace=0.12)
    thumb = fig.add_subplot(gs[0, 0])
    gloss = fig.add_subplot(gs[1, 0]); gloss.axis("off")
    ax = fig.add_subplot(gs[:, 1])

    # thumbnail
    thumb.plot(pf, Ef, color=CE, lw=2.2); thumb.plot(pf, Lf, color=CL, lw=2.2)
    thumb.axhline(Eh, color="#888888", lw=1, ls="--")
    y0, y1 = Eh + yL[p >= lo].min() - 0.01, Eh + yE.max() + 0.01
    thumb.add_patch(Rectangle((lo, y0), hi - lo, y1 - y0, fill=False, ec="#2b6cb0", lw=2.5))
    thumb.annotate("the window on the right", (hi, y1), xytext=(40, 40), textcoords="offset points",
                   color="#2b6cb0", arrowprops=dict(arrowstyle="->", color="#2b6cb0", lw=2))
    thumb.set_xlim(0, 1); thumb.set_ylim(Eh - 0.25*(n - Eh), n + 0.05*(n - Eh))
    thumb.set_title(f"n = {n}, all of [0, 1]:  E (black) and L (red)", loc="left")
    thumb.set_xlabel("p")
    gloss.text(0.0, 1.0, GLOSSARY, va="top", ha="left", fontsize=16.5, linespacing=1.45,
               transform=gloss.transAxes)

    # the window: certificates, slack, E, L
    for cv, col in CC.items():
        k = int(np.argmin(np.abs(c - cv)))
        ax.plot(p, yW[k], color=col, lw=2.0, alpha=0.75, zorder=2)
    ax.fill_between(p, yL, yE, color=CS, alpha=0.75, lw=0, zorder=3)
    ax.plot(p, yE, color=CE, lw=3.4, zorder=5)
    ax.plot(p, yL, color=CL, lw=3.4, zorder=6)
    ax.axhline(0, color="#888888", lw=1.4, ls="--", zorder=1)

    tw = ties[(ties > 0.5) & (ties >= lo) & (ties <= hi)]
    ax.plot(tw, [at(yE, v) for v in tw], "|", color=CE, ms=22, mew=2.2, zorder=7)
    cw = cp[(cp >= lo) & (cp <= hi)]
    ax.plot(cw, [at(yE, v) for v in cw], "o", ms=15, mfc="white", mec=CE, mew=2.8, zorder=8)
    sww = sw[(sw >= lo) & (sw <= hi)]
    for k, s in enumerate(sww, 1):
        ax.axvline(s, color=CL, lw=1.3, ls=":", zorder=0)
        ax.plot([s], [at(yL, s)], "D", ms=12, color=CL, zorder=8)
    ax.plot([0.5], [0], "s", ms=13, color="#444444", zorder=9)

    ymin, ymax = -0.046, 0.092
    ax.set_xlim(lo, hi); ax.set_ylim(ymin, ymax)
    ax.set_xlabel("p"); ax.set_ylabel("value $-$ E(n, 1/2)")
    sec = ax.secondary_xaxis("top", functions=(lambda v: n*(v - 0.5), lambda x: 0.5 + x/n))
    sec.set_xlabel("x = n(p $-$ 1/2)")
    ax.set_title(f"n = {n}, close up above p = 1/2", loc="left", pad=60)

    # ---- labels: short names on the plot, definitions in the glossary (positions tuned for n = 12)
    B = dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.85)
    def lab(txt, xy, xytext, col, ha="center", va="center"):
        ax.annotate(txt, xy, xytext=xytext, color=col, ha=ha, va=va, fontsize=20, zorder=20, bbox=B,
                    arrowprops=None if xy is None else dict(arrowstyle="->", color=col, lw=2, shrinkB=7))
    def name(txt, x, y, col, ha):
        ax.text(x, y, txt, color=col, ha=ha, va="center", fontsize=21, zorder=20, bbox=B)
    pair = lambda v: pairs.get(round(float(v), 12), ("?", "?"))
    sw1, sw2 = sww[0], sww[1]

    # the certificates, named on stretches where they are clear of L
    name("W$_{5.75}$", 0.5095, -0.024, CC[5.75], "left")
    name("W$_{6.25}$", 0.5505, -0.030, CC[6.25], "right")
    name("W$_{6.75}$", 0.5305, -0.018, CC[6.75], "right")
    name("W$_{7.25}$", 0.5650, -0.030, CC[7.25], "right")
    # E, L, slack
    pe = 0.5878
    lab("E(n, p)", (pe, at(yE, pe)), (0.5765, 0.069), CE, ha="right")
    lab("L(p) = max$_c$ W$_c$(p)", (pe, at(yL, pe)), (0.5915, -0.022), CL, ha="right")
    lab("slack E $-$ L", (0.5852, 0.5*(at(yE, 0.5852) + at(yL, 0.5852))), (0.5835, 0.085), "#555555")
    lab("slack", (0.5428, 0.5*(at(yE, 0.5428) + at(yL, 0.5428))), (0.5480, 0.008), "#555555")
    # p = 1/2
    lab("p = 1/2:  E = L", (0.5, 0.0), (0.4975, -0.040), "#444444", ha="left")
    # active certificates
    lab("active: W$_{6.25}$   (E = L, slack 0)", (0.5135, at(yL, 0.5135)), (0.4985, 0.045), CC[6.25], ha="left")
    lab("active: W$_{6.75}$", (0.5615, at(yL, 0.5615)), (0.5555, 0.027), CC[6.75])
    # switch points
    lab("switch point p$_1$ (HALF)", (sw1, at(yL, sw1)), (0.5385, -0.006), CL, ha="right")
    lab("switch point p$_2$ (INT)", (sw2, at(yL, sw2)), (0.5740, -0.010), CL, ha="right")
    # tie points and cusps
    c1 = cw[np.argmin(np.abs(cw - 0.539))]
    lab("cusp (%d,%d)" % pair(c1), (c1, at(yE, c1)), (0.5315, 0.040), CE)
    two = sorted(cw[(cw > 0.57) & (cw < 0.585)])
    lab("cusps (%d,%d) and (%d,%d)" % (pair(two[0]) + pair(two[1])), (two[0], at(yE, two[0])),
        (0.5700, 0.062), CE, ha="right")
    lab("", (two[1], at(yE, two[1])), (0.5700, 0.062), CE)
    t1 = tw[np.argmin(np.abs(tw - 0.5562))]
    lab("tie point (%d,%d)" % pair(t1), (t1, at(yE, t1)), (0.5470, 0.052), CE)
    # margin of L
    pm = 0.5672
    ax.annotate("", (pm, at(yL, pm)), xytext=(pm, 0), zorder=19,
                arrowprops=dict(arrowstyle="<->", color="#2b6cb0", lw=2.6, shrinkA=0, shrinkB=0))
    name("margin of L", pm + 0.0012, 0.022, "#2b6cb0", "left")

    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, f"L_glossary_n{n:03d}.png")
    fig.savefig(path, dpi=100, bbox_inches="tight")
    print(path)

if __name__ == "__main__":
    main()
