"""plotting/simplex/triangle_n2.py -- the n = 2 simplex: chambers, the sorted sliver, E's level sets,
B = {E <= E(1/2)}, {L <= L(1/2)}, the certificates, and the binomial curve.

    .venv/bin/python plotting/simplex/triangle_n2.py [--out DIR]

x = (f0, f1, f2) = (q^2, 2pq, p^2) runs from e0 (p = 0) to e2 (p = 1).
  chambers   the six regions with a fixed ordering of (f0, f1, f2), cut by the walls f_i = f_j.  The
             curve crosses a wall at each tie point.  The two chambers where f1 is the SMALLEST are
             never entered (binomial masses are unimodal) and no certificate ranks f1 lowest: hatched.
  sliver     the chamber f0 <= f1 <= f2 where the SORTED vector lives.  sort(f_p) folds the curve
             into it, bouncing off its two walls at the tie points (dotted).
  E          linear on each chamber, so its level sets are hexagons around the centre (thin grey);
             B = {E <= E(1/2) = 5/4} is the conjecture's forbidden body: the curve touches it only at
             p = 1/2, at its vertex (1/4, 1/2, 1/4).
  L          {L <= L(1/2)} is cut out by one half-plane per certificate (W_-1/4 and W_9/4 are dominated
             by W_1/4 and W_7/4 everywhere and are not drawn) W_c(x) = <a_c, x> <= 5/4.
             Since L <= E it contains B.  Near the touch point its two edges ARE B's two facets
             (W_3/4 and W_5/4, the rankings 1>0>2 and 1>2>0); away from it the certificates are weaker
             than B's facets, and below the apex it opens down to the edge f1 = 0.  (★) is the
             statement that the curve avoids this larger region too.
Output: simplex_figures/triangle_n2.png.
"""
import argparse, math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _simplex import (OUT, TRI, xy, chamber, E_sublevel, L_sublevel, certificates, binomial,
                      line_in_simplex)

SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
C_B, C_L, C_SLIV, C_HATCH = "#3a7dc9", "#d1495b", "#f6e7b8", "#b8bec6"
FRAC = {-0.25: "−1/4", 0.25: "1/4", 0.75: "3/4", 1.25: "5/4", 1.75: "7/4", 2.25: "9/4"}

def poly_patch(poly, **kw):
    return Polygon(xy(np.array(poly)), closed=True, **kw)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    n, E0 = 2, 1.25

    plt.rcParams.update({"font.size": 18})
    fig, ax = plt.subplots(figsize=(22, 18), dpi=100)
    ax.set_aspect("equal"); ax.axis("off")

    # chambers: hatch the two where f1 is smallest; tint the sorted sliver
    for order in ((0, 2, 1), (2, 0, 1)):
        ax.add_patch(poly_patch(chamber(order), fc="none", ec=C_HATCH, hatch="//", lw=0, zorder=1))
    ax.add_patch(poly_patch(chamber((2, 1, 0)), fc=C_SLIV, ec="none", zorder=1))

    # E level sets (hexagons), then B and {L <= L(1/2)}
    for t in np.arange(1.05, 2.0, 0.1):
        if abs(t - E0) > 1e-9:
            ax.add_patch(poly_patch(E_sublevel(n, t), fc="none", ec="#9aa0a6", lw=1.0, zorder=2))
    ax.add_patch(poly_patch(L_sublevel(n, E0), fc=C_L, alpha=0.10, ec="none", zorder=3))
    ax.add_patch(poly_patch(L_sublevel(n, E0), fc="none", ec=C_L, lw=3.0, ls=(0, (6, 3)), zorder=6))
    ax.add_patch(poly_patch(E_sublevel(n, E0), fc=C_B, alpha=0.28, ec=C_B, lw=3.0, zorder=4))

    # walls f_i = f_j: from the vertex e_k through the centre to the opposite edge's midpoint
    walls = {(0, 1): 2, (0, 2): 1, (1, 2): 0}
    for (i, j), k in walls.items():
        mid = np.zeros(3); mid[i] = mid[j] = 0.5
        P = xy(np.array([np.eye(3)[k], mid]))
        ax.plot(P[:, 0], P[:, 1], color="#555555", lw=1.4, ls=(0, (2, 3)), zorder=5)
        if k != 1: continue                  # the other two walls are named by their tie points
        lab = xy(mid) + (xy(mid) - xy(np.eye(3)[k]))*0.075
        ax.text(*lab, f"f{i} = f{j}".translate(SUB), ha="center", va="center", fontsize=17,
                color="#444444", rotation=0)

    # certificate lines W_c = 5/4, labelled at their right/lower end inside the triangle
    certs = certificates(n)
    ccol = {-0.25: "#9b8fc7", 0.25: "#7b6fb0", 0.75: "#2a9d8f", 1.25: "#e9a03b", 1.75: "#3a7dc9", 2.25: "#6fa8dc"}
    for c, av in certs.items():
        if c in (-0.25, 2.25): continue      # dominated: W_1/4 - W_-1/4 = (0,1,1) >= 0, never active
        seg = line_in_simplex(n, av, E0)
        if len(seg) < 2: continue
        P = xy(np.array(seg))
        main_pair = c in (0.75, 1.25)
        ax.plot(P[:, 0], P[:, 1], color=ccol[c], lw=2.4 if main_pair else 1.5,
                alpha=0.95 if main_pair else 0.7, zorder=7)
        end = P[np.argmax(P[:, 1])]
        ax.annotate(f"W$_{{{FRAC[c]}}}$ = 5/4", end, xytext=(10 if end[0] > 0.5 else -10, 4),
                    textcoords="offset points", ha="left" if end[0] > 0.5 else "right",
                    color=ccol[c], fontsize=17, zorder=12)

    # the binomial curve, coloured by p, and its sorted (folded) copy in the sliver
    p = np.linspace(0, 1, 1201)
    f = binomial(n, p)
    P = xy(f)
    segs = np.stack([P[:-1], P[1:]], axis=1)
    cmap = plt.get_cmap("coolwarm")
    lc = LineCollection(segs, cmap=cmap, lw=5.5, zorder=9, capstyle="round")
    lc.set_array(0.5*(p[:-1] + p[1:])); lc.set_clim(0, 1)
    ax.add_collection(lc)
    S = xy(np.sort(f, axis=1))
    ax.plot(S[:, 0], S[:, 1], color="#333333", lw=2.2, ls=(0, (1, 2)), zorder=8)
    for pt in np.arange(0.1, 1.0, 0.1):
        if abs(pt - 0.5) < 1e-9: continue
        q = xy(binomial(n, [pt]))[0]
        ax.plot(*q, "o", ms=6, color="black", zorder=10)
        off = (-14, 10) if pt < 0.5 else (14, 10)
        ax.annotate(f"{pt:.1f}", q, xytext=off, textcoords="offset points", fontsize=15,
                    ha="right" if pt < 0.5 else "left", color="#333333", zorder=12)

    # tie points: the crossings of the walls
    ties = {(0, 1): 1/3, (0, 2): 1/2, (1, 2): 2/3}
    for (i, j), pt in ties.items():
        q = xy(binomial(n, [pt]))[0]
        ax.plot(*q, "D", ms=13, mfc="white", mec="black", mew=2.4, zorder=11)
    ax.annotate("tie point (0,1), p = 1/3,\non the wall f₀ = f₁", xy(binomial(n, [1/3]))[0], xytext=(-150, 70),
                textcoords="offset points", ha="right", fontsize=17, zorder=12,
                arrowprops=dict(arrowstyle="->", lw=1.6))
    ax.annotate("tie point (1,2), p = 2/3,\non the wall f₁ = f₂", xy(binomial(n, [2/3]))[0], xytext=(150, 70),
                textcoords="offset points", ha="left", fontsize=17, zorder=12,
                arrowprops=dict(arrowstyle="->", lw=1.6))
    touch = xy(np.array([0.25, 0.5, 0.25]))
    ax.annotate("p = 1/2, tie point (0,2), on the wall f₀ = f₂:\nthe curve touches B and {L ≤ L(1/2)}\nat their shared vertex (1/4, 1/2, 1/4)",
                touch, xytext=(-0.13, 0.80), textcoords="data", ha="left", va="center", fontsize=17, zorder=12,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.9),
                arrowprops=dict(arrowstyle="->", lw=1.8))

    # regions named
    def tag(txt, bary, col, fontsize=18):
        ax.text(*xy(np.array(bary)), txt, color=col, ha="center", va="center", fontsize=fontsize, zorder=12,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.85))
    tag("B = {E ≤ E(1/2)}", [0.29, 0.42, 0.29], C_B)
    tag("{L ≤ L(1/2)}  ⊇  B", [0.465, 0.07, 0.465], C_L)
    note = dict(textcoords="data", fontsize=17, ha="center", va="center", zorder=12,
                arrowprops=dict(arrowstyle="->", lw=1.6, color="#6b7280"))
    hx = xy(np.array([0.5, -0.15, 0.5]))
    ax.annotate("f₁ smallest (hatched): never visited by any f$_p$,\nand no certificate ranks f₁ last",
                xy(np.array([0.84, 0.03, 0.13])), xytext=hx, color="#6b7280", **note)
    ax.annotate("", xy(np.array([0.13, 0.03, 0.84])), xytext=hx + np.array([0.13, 0.025]), **note)
    ax.annotate("sorted sliver\nf₀ ≤ f₁ ≤ f₂\n(dotted: sort(f$_p$),\nthe ordered version)",
                xy(np.array([0.17, 0.30, 0.53])), xytext=(1.03, 0.36),
                color="#8a6d1c", ha="left", va="center", fontsize=17, zorder=12,
                arrowprops=dict(arrowstyle="->", lw=1.6, color="#8a6d1c"))

    # triangle and vertices
    ax.add_patch(Polygon(TRI, closed=True, fc="none", ec="black", lw=2.5, zorder=13))
    vt = {0: ("e₀ = (1, 0, 0)\np = 0", (-20, -10), "right"), 1: ("e₁ = (0, 1, 0)", (0, 22), "center"),
          2: ("e₂ = (0, 0, 1)\np = 1", (20, -10), "left")}
    for k, (txt, off, ha) in vt.items():
        ax.annotate(txt, TRI[k], xytext=off, textcoords="offset points", ha=ha, va="center", fontsize=19)
    ax.plot(*xy(np.array([1/3, 1/3, 1/3])), "+", ms=16, mew=2, color="#555555", zorder=12)

    # key
    hand = [Line2D([], [], color=cmap(0.85), lw=5, label="binomial curve f$_p$ = (q², 2pq, p²), coloured by p (blue 0, red 1)"),
            Line2D([], [], color="#9aa0a6", lw=1.2, label="level sets of E (hexagons; E is linear on each chamber)"),
            Line2D([], [], color="#555555", lw=1.4, ls=(0, (2, 3)), label="walls f$_i$ = f$_j$ (tie points lie on them)"),
            Line2D([], [], color="#2a9d8f", lw=2.4, label="certificate lines W$_c$(x) = 5/4  (W$_c$ = ⟨a$_c$, x⟩, a$_c$[k] = n + 1/2 − 2|k − c|)"),
            Line2D([], [], color="#333333", lw=2.2, ls=(0, (1, 2)), label="sort(f$_p$): the ordered version, folded into the sliver"),
            Line2D([], [], ls="", marker="D", ms=11, mfc="white", mec="black", mew=2.2, label="tie points: the curve crosses a wall")]
    ax.legend(handles=hand, loc="upper center", fontsize=16, frameon=False, bbox_to_anchor=(0.5, 0.0), ncol=2)
    ax.set_title("n = 2: the simplex of mass vectors (f₀, f₁, f₂)", loc="left", fontsize=24)
    ax.set_xlim(-0.14, 1.24); ax.set_ylim(-0.19, math.sqrt(3)/2 + 0.06)

    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, "triangle_n2.png")
    fig.savefig(path, dpi=100, bbox_inches="tight")
    print(path)

if __name__ == "__main__":
    main()
