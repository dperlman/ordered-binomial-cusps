"""plotting/simplex/kaleidoscope_n2.py -- E's and L's folded paths in the n = 2 triangle.

    .venv/bin/python plotting/simplex/kaleidoscope_n2.py [--out DIR]

Every rearrangement of f_p is a mirror image of the binomial curve: six copies, one per chamber (the
kaleidoscope).  E and L each pick one copy at every p:
  E's path    y_E(p) = sort(f_p): the copy in the sliver y0 <= y1 <= y2.  Continuous; it BOUNCES off a
              wall at each tie point.
  L's path    y_L(p) = f_p rearranged by the ordering of L's active certificate (rank by distance from
              its centre c; the mass ranked k-th lowest goes to slot k).  It JUMPS at L's switch points.
The height h(y) = <(0, 1, 2), y> = y1 + 2 y2 (= 2 x horizontal position here) gives E(p) = h(y_E).
L scores certificate c with the formula a_c(k) = n + 1/2 - 2|k - c|, not with the plain ranks sigma_c:
the formula ranks by distance on the whole integer line, so positions outside 0..n ("phantoms") use up
rank slots, and W_c = <a_c, f_p> = <sigma_c, f_p> - Delta_c with the PHANTOM DEFICIT
Delta_c = E_p[(2c - n - 1/2 - K)^+] >= 0 (c > n/2).  The central certificates (c = n/2 +- 1/4) have none.
At n = 2 every ordering a binomial vector ever takes is one of the distance orderings, so with plain
ranks the best certificate would always be the true ordering (L* = E).  All of n = 2's slack comes from
the phantom deficit:
  1/2 < p < 2/3        L uses W_5/4 (ordering 1 > 2 > 0, true, no phantoms): paths coincide, E = L.
  2/3 < p < 1/sqrt(2)  the true ordering is now 2 > 1 > 0, W_7/4's -- but W_7/4 is scored low by its
                       deficit q^2, so L keeps W_5/4.  L's path runs on the mirror image of E's path,
                       and E - L is the height gap between the two paths.
  p = 1/sqrt(2)        W_7/4's formula score catches up: L switches and its path jumps back.
  p > 1/sqrt(2)        the paths coincide again, and L = E - q^2: the deficit itself, not a point.
Since sort(f_{1-p}) = sort(f_p) and y_L(1-p) = y_L(p), both paths are traced out and back; only
p in [1/2, 1] is drawn.  Output: simplex_figures/kaleidoscope_n2.png.
"""
import argparse, itertools, math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon, Rectangle
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _simplex import OUT, TRI, xy, chamber, binomial, certificates

N = 2
CE, CL, CSL, CCOPY = "#111111", "#d1495b", "#f6e7b8", "#b8bec6"
SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")

def ranking(c, n=N):
    """sigma_c: rank of each k by distance from c (nearest = n)."""
    order = np.argsort(np.abs(np.arange(n + 1) - c))
    r = np.empty(n + 1, int); r[order] = n - np.arange(n + 1)
    return r

def paths(p):
    f = binomial(N, p)
    certs = {c: a for c, a in certificates(N).items() if c not in (-0.25, 2.25)}   # the dominated two never win
    cs = np.array(list(certs)); Avec = np.array(list(certs.values()))
    act = cs[np.argmax(f @ Avec.T, axis=1)]
    yL = np.empty_like(f)
    for t, c in enumerate(act):
        yL[t, ranking(c)] = f[t]
    yE = np.sort(f, axis=1)
    h = np.arange(N + 1)
    return dict(f=f, act=act, yE=yE, yL=yL, E=yE @ h, L=np.max(f @ Avec.T, axis=1),
                W54=f @ certs[1.25], W74=f @ certs[1.75])

def draw_triangle(ax, P, p, zoom):
    """Kaleidoscope, sliver, walls, height lines, and the two paths."""
    ax.add_patch(Polygon(xy(np.array(chamber((2, 1, 0)))), closed=True, fc=CSL, ec="none", zorder=0))
    # height level lines: in this embedding the horizontal coordinate is y1/2 + y2 = h/2, so the
    # level sets of h are vertical lines and the height is twice the horizontal position
    for hv in np.arange(0.1, 2.0, 0.1):
        xv = hv/2
        ax.plot([xv, xv], [0, math.sqrt(3)*min(xv, 1 - xv)], color="#d7dbe0", lw=1.0, zorder=1)
    # walls
    for (i, j), k in {(0, 1): 2, (0, 2): 1, (1, 2): 0}.items():
        mid = np.zeros(3); mid[i] = mid[j] = 0.5
        Q = xy(np.array([np.eye(3)[k], mid]))
        ax.plot(Q[:, 0], Q[:, 1], color="#555555", lw=1.3, ls=(0, (2, 3)), zorder=2)
    # the six mirror copies of the curve
    pf = np.linspace(0, 1, 1201); ff = binomial(N, pf)
    for perm in itertools.permutations(range(3)):
        C = xy(ff[:, perm])
        ident = perm == (0, 1, 2)
        ax.plot(C[:, 0], C[:, 1], color="#7d8590" if ident else CCOPY, lw=1.8 if ident else 1.4,
                ls="-" if ident else "-", zorder=3, alpha=1.0)
    # E's path
    YE = xy(P["yE"]); ax.plot(YE[:, 0], YE[:, 1], color=CE, lw=5.0 if zoom else 3.5, zorder=6, solid_capstyle="round")
    # L's path, split at the switch points; jumps dashed
    YL = xy(P["yL"]); act = P["act"]
    cuts = np.flatnonzero(np.diff(act) != 0) + 1
    for seg in np.split(np.arange(len(p)), cuts):
        ax.plot(YL[seg, 0], YL[seg, 1], color=CL, lw=2.4 if zoom else 1.8, ls=(0, (5, 3)), zorder=7)
    for k in cuts:
        ax.annotate("", YL[k], xytext=YL[k - 1], zorder=8,
                    arrowprops=dict(arrowstyle="->", color=CL, lw=2.2, ls="-", shrinkA=0, shrinkB=0))
    ax.add_patch(Polygon(TRI, closed=True, fc="none", ec="black", lw=2.2, zorder=9))
    ax.set_aspect("equal"); ax.axis("off")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    ptie, psw = 2/3, 1/math.sqrt(2)
    p = np.unique(np.concatenate([np.linspace(0.5, 1, 6001), [ptie, psw - 1e-12, psw + 1e-12]]))
    P = paths(p)

    plt.rcParams.update({"font.size": 17})
    fig = plt.figure(figsize=(32, 22), dpi=100)
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.15], height_ratios=[1.25, 1], hspace=0.12, wspace=0.05)

    # (a) the whole kaleidoscope
    ax = fig.add_subplot(gs[0, 0]); draw_triangle(ax, P, p, False)
    for k, (txt, off, ha) in {0: ("e₀", (-14, -8), "right"), 1: ("e₁", (0, 14), "center"), 2: ("e₂", (14, -8), "left")}.items():
        ax.annotate(txt, TRI[k], xytext=off, textcoords="offset points", ha=ha, fontsize=20)
    ax.annotate("f$_p$ itself", xy(binomial(N, [0.2]))[0], xytext=(-90, 30), textcoords="offset points",
                color="#5f6670", fontsize=16, arrowprops=dict(arrowstyle="->", color="#5f6670"))
    ax.text(*xy(np.array([0.08, 0.30, 0.62])), "sliver", color="#8a6d1c", fontsize=17, ha="center")
    zx0, zx1, zy0, zy1 = 0.575, 0.875, 0.185, 0.465
    ax.add_patch(Rectangle((zx0, zy0), zx1 - zx0, zy1 - zy0, fill=False, ec="#2b6cb0", lw=2.2, zorder=10))
    ax.set_title("(a) the six mirror copies of the curve (grey): every rearrangement of f$_p$", loc="left", fontsize=20)

    # (b) close up: the late bounce and the jump
    bx = fig.add_subplot(gs[0, 1]); draw_triangle(bx, P, p, True)
    bx.set_xlim(zx0, zx1); bx.set_ylim(zy0, zy1); bx.axis("on")
    hs = np.round(np.arange(1.2, 1.75, 0.1), 1)
    bx.set_xticks(hs/2); bx.set_xticklabels([f"{v:g}" for v in hs]); bx.set_yticks([])
    bx.xaxis.tick_top(); bx.xaxis.set_label_position("top")
    bx.set_xlabel("height h = y₁ + 2y₂ = 2 × horizontal position (grey vertical lines)", fontsize=16)
    for s in bx.spines.values(): s.set_color("#2b6cb0"); s.set_linewidth(2)
    def at(pv, key):
        return xy(P[key][int(np.argmin(np.abs(p - pv)))])
    A = dict(textcoords="offset points", fontsize=17, zorder=20,
             bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.9))
    bx.plot(*at(0.5, "yE"), "s", ms=12, color="black", zorder=12)
    bx.annotate("p = 1/2: both paths start here,\nat (1/4, 1/4, 1/2), height 5/4", at(0.5, "yE"), xytext=(-40, 20), ha="right", **A)
    bx.plot(*at(ptie, "yE"), "D", ms=13, mfc="white", mec="black", mew=2.4, zorder=12)
    bx.annotate("tie point (1,2), p = 2/3:\nE's path bounces off the wall y₁ = y₂", at(ptie, "yE"),
                xytext=(40, 40), ha="left", arrowprops=dict(arrowstyle="->", lw=1.6), **A)
    k0 = int(np.argmin(np.abs(p - (psw - 1e-12)))); k1 = k0 + 1
    bx.plot(*xy(P["yL"][k0]), "o", ms=11, color=CL, zorder=12)
    bx.plot(*xy(P["yL"][k1]), "o", ms=11, mfc="white", mec=CL, mew=2.4, zorder=12)
    bx.annotate("switch point p = 1/√2 ≈ 0.707: W$_{7/4}$'s formula\nscore catches up with W$_{5/4}$'s, and L's path\njumps back onto E's",
                xy(P["yL"][k0]), xytext=(70, 10), ha="left", color=CL, arrowprops=dict(arrowstyle="->", lw=1.6, color=CL), **A)
    mid = 0.5*(ptie + psw)
    bx.annotate("2/3 < p < 0.707: the true ordering is now 2 > 1 > 0,\nwhich is W$_{7/4}$'s, but the formula scores W$_{7/4}$\nlow by its phantom deficit q², so L keeps\nW$_{5/4}$'s ordering 1 > 2 > 0: its path runs on\nthe MIRROR IMAGE of E's path across the wall",
                at(mid, "yL"), xytext=(-40, 10), ha="right", color=CL, arrowprops=dict(arrowstyle="->", lw=1.6, color=CL), **A)
    bx.annotate("1/2 < p < 2/3: L's ordering is the true one,\nthe paths coincide, E = L",
                at(0.6, "yE"), xytext=(-60, 10), ha="right", arrowprops=dict(arrowstyle="->", lw=1.6), **A)
    bx.annotate("p > 0.707: the paths coincide again,\nheading for e₂ (p = 1); L = E − q²,\nthe phantom deficit (panel c)",
                at(0.74, "yE"), xytext=(-40, -110), ha="right", arrowprops=dict(arrowstyle="->", lw=1.6), **A)
    bx.set_title("(b) close up (blue box in a): L switches late, so its path crosses the wall y₁ = y₂", loc="left", fontsize=20, pad=48)

    # (c) heights, and the two certificate scores competing above 1/2
    cx = fig.add_subplot(gs[1, :])
    C54, C74 = "#2a9d8f", "#e9a03b"
    cx.plot(p, P["W54"], color=C54, lw=1.8, label="W$_{5/4}$ = ⟨a$_{5/4}$, f$_p$⟩: ordering 1 > 2 > 0, no phantom deficit")
    cx.plot(p, P["W74"], color=C74, lw=1.8, label="W$_{7/4}$ = ⟨a$_{7/4}$, f$_p$⟩: ordering 2 > 1 > 0 (the true one for p > 2/3), scored q² low")
    cx.plot(p, P["E"], color=CE, lw=3.4, label="E(p) = h(y$_E$), the height of E's path")
    cx.plot(p, P["L"], color=CL, lw=3.4, ls=(0, (5, 2)), label="L(p) = max$_c$ W$_c$: the higher of the two")
    w = (p > ptie) & (p < psw)
    cx.fill_between(p, P["L"], P["E"], where=w, color="#9aa5b1", alpha=0.6, lw=0,
                    label="L switches late: its path is on the mirror copy, E − L = the height gap between the paths")
    cx.fill_between(p, P["L"], P["E"], where=p > psw, color=CL, alpha=0.15, lw=0,
                    label="paths coincide; E − L = q², W$_{7/4}$'s phantom deficit (not a point in the triangle)")
    cx.axvline(ptie, color="black", lw=1, ls=":"); cx.axvline(psw, color=CL, lw=1, ls=":")
    cx.text(ptie, 1.235, "tie point 2/3:  \nE's ordering changes  ", va="bottom", ha="right", fontsize=15)
    cx.text(psw, 1.235, " switch point 1/√2:\n W$_{7/4}$ overtakes W$_{5/4}$", va="bottom", fontsize=15, color=CL)
    cx.set_xlim(0.5, 1); cx.set_ylim(1.2, 2.02)
    cx.set_xlabel("p   (p < 1/2 retraces the same paths: sort(f$_{1-p}$) = sort(f$_p$), y$_L$(1−p) = y$_L$(p))")
    cx.set_ylabel("height h = y₁ + 2 y₂")
    cx.legend(loc="upper left", fontsize=15.5, frameon=False)
    cx.set_title("(c) as heights against p: L is the higher of the two certificate scores; at n = 2 all of E − L "
                 "comes from W$_{7/4}$'s phantom deficit", loc="left", fontsize=20)

    hand = [Line2D([], [], color=CE, lw=4, label="E's path y$_E$(p) = sort(f$_p$)"),
            Line2D([], [], color=CL, lw=2.4, ls=(0, (5, 3)), label="L's path y$_L$(p): f$_p$ rearranged by L's active certificate's ordering"),
            Line2D([], [], color=CCOPY, lw=1.6, label="mirror copies of the curve (rearrangements of f$_p$)"),
            Line2D([], [], color="#d7dbe0", lw=1.4, label="height level lines h = y₁ + 2 y₂ = const"),
            Line2D([], [], color="#555555", lw=1.3, ls=(0, (2, 3)), label="walls y$_i$ = y$_j$")]
    fig.legend(handles=hand, loc="upper center", ncol=3, fontsize=17, frameon=False, bbox_to_anchor=(0.5, 0.965))
    fig.suptitle("n = 2: E and L as two folded paths through the same kaleidoscope.\nEvery ordering the binomial takes "
                 "is in L's menu; L strays only because its formula marks one ordering down (the phantom deficit).",
                 x=0.01, y=1.01, ha="left", fontsize=24)
    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, "kaleidoscope_n2.png")
    fig.savefig(path, dpi=100, bbox_inches="tight")
    print(path)

if __name__ == "__main__":
    main()
