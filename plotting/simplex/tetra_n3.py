"""plotting/simplex/tetra_n3.py -- the n = 3 simplex (a tetrahedron): chambers, the sorted sliver,
B = {E <= E(1/2)}, {L <= L(1/2)}, the certificates and the binomial curve.  The n = 2 triangle, one
dimension up.

    .venv/bin/python plotting/simplex/tetra_n3.py [--out DIR]

Writes simplex_figures/tetra_n3.html (interactive, plotly: drag to rotate, click legend entries to
show/hide layers; some start hidden) and simplex_figures/tetra_n3.png (three fixed views, matplotlib).

x = (q^3, 3pq^2, 3p^2q, p^3) runs from e0 (p = 0) to e3 (p = 1).  E(1/2) = L(1/2) = 2.
  B          {E <= 2}: one facet per permutation -- 14 vertices, 24 triangles (a tetrakis hexahedron,
             the polar of the permutohedron).  The curve touches it only at p = 1/2, where TWO pairs tie
             at once ((0,3) and (1,2)), so the touch point lies on the line where two walls meet.
  {L <= 2}   one half-space per certificate W_c (c = -1/4 .. 13/4); contains B, reaches the boundary
             of the tetrahedron, and the curve stays outside it -- (★) for n = 3.
  sliver     the chamber f0 <= f1 <= f2 <= f3 (1/24 of the volume) where sort(f_p) lives.
  chambers   24 orderings; only the 8 unimodal ones can hold a binomial vector.  The other 16 are a
             hidden layer in the HTML.
  ties       above 1/2: (1,3) at p = sqrt(3)/(1+sqrt(3)) = 0.634 -- the one cusp of n = 3 -- and (2,3) at
             p = 3/4; their mirrors (0,2) and (0,1) below.
"""
import argparse, itertools, math, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _simplex import OUT, binomial, certificates, perms
from _simplex3d import TET, xyz, polytope, face_on_plane, edges, triangles

N, E0 = 3, 2.0
CEN = np.full(4, 0.25)
C_B, C_L, C_SLIV = "#3a7dc9", "#d1495b", "#d9a400"
VNAME = ["e₀ = (1,0,0,0)  p = 0", "e₁ = (0,1,0,0)", "e₂ = (0,0,1,0)", "e₃ = (0,0,0,1)  p = 1"]
TIES = {(0, 1): 0.25, (0, 2): 1/(1 + math.sqrt(3)), (0, 3): 0.5, (1, 2): 0.5,
        (1, 3): math.sqrt(3)/(1 + math.sqrt(3)), (2, 3): 0.75}
CUSPS = {(0, 2), (1, 3)}

def chamber(order):
    """Polytope of x[order[0]] >= x[order[1]] >= ... ."""
    hs = []
    for hi, lo in zip(order[:-1], order[1:]):
        a = np.zeros(4); a[lo] = 1; a[hi] = -1
        hs.append((a, 0.0))
    inner = np.zeros(4); inner[list(order)] = [0.4, 0.3, 0.2, 0.1]
    return polytope(hs, inner)

def unimodal(order):
    r = np.empty(4); r[list(order)] = [3, 2, 1, 0]          # rank of each k
    return not any(r[k] < r[k - 1] and r[k] < r[k + 1] for k in (1, 2))

def wall(i, j):
    """The wall x_i = x_j inside the tetrahedron: triangle through e_k, e_l and mid(e_i, e_j)."""
    k, l = [m for m in range(4) if m not in (i, j)]
    mid = np.zeros(4); mid[i] = mid[j] = 0.5
    return xyz(np.array([np.eye(4)[k], np.eye(4)[l], mid]))

def geometry():
    g = {}
    g["B"] = polytope([(s, E0) for s in perms(N)], CEN)
    certs = certificates(N)
    g["L"] = polytope([(a, E0) for a in certs.values()], CEN)
    g["sliver"] = chamber((3, 2, 1, 0))
    g["never"] = [chamber(o) for o in itertools.permutations(range(4)) if not unimodal(o)]
    g["levels"] = {t: polytope([(s, t) for s in perms(N)], CEN) for t in (2.25, 2.5, 2.75)}
    g["certs"] = {}
    for c, a in certs.items():                          # the plane W_c = 2 cut to the tetrahedron
        V, F = polytope([(a, E0)], CEN) if a @ CEN < E0 else polytope([(-a, -E0)], CEN)
        g["certs"][c] = face_on_plane(V, F, a, E0)
    p = np.linspace(0, 1, 1201)
    f = binomial(N, p)
    g["p"], g["curve"], g["sorted"] = p, xyz(f), xyz(np.sort(f, axis=1))
    return g

# ------------------------------------------------------------------------------------------- plotly
def html(g, path):
    import plotly.graph_objects as go
    tr = []
    def segs(seglist):
        xs, ys, zs = [], [], []
        for s in seglist:
            xs += [s[0][0], s[1][0], None]; ys += [s[0][1], s[1][1], None]; zs += [s[0][2], s[1][2], None]
        return xs, ys, zs
    def mesh(V, F, color, opacity, name, group, show=True, legend=True):
        T = np.array(triangles(F))
        tr.append(go.Mesh3d(x=V[:, 0], y=V[:, 1], z=V[:, 2], i=T[:, 0], j=T[:, 1], k=T[:, 2], color=color,
                            opacity=opacity, flatshading=True, name=name, legendgroup=group,
                            showlegend=legend, visible=True if show else "legendonly", hoverinfo="name"))
    def wire(seglist, color, width, name, group, show=True, legend=False, dash="solid"):
        x, y, z = segs(seglist)
        tr.append(go.Scatter3d(x=x, y=y, z=z, mode="lines", line=dict(color=color, width=width, dash=dash),
                               name=name, legendgroup=group, showlegend=legend,
                               visible=True if show else "legendonly", hoverinfo="skip"))

    # tetrahedron
    wire([TET[[i, j]] for i, j in itertools.combinations(range(4), 2)], "black", 4, "tetrahedron", "tet", legend=True)
    tr.append(go.Scatter3d(x=TET[:, 0]*1.12, y=TET[:, 1]*1.12, z=TET[:, 2]*1.12, mode="text", text=VNAME,
                           textfont=dict(size=15), showlegend=False, hoverinfo="skip", legendgroup="tet"))
    # B and {L <= 2}
    V, F = g["B"]; mesh(V, F, C_B, 0.45, "B = {E ≤ E(1/2)}", "B"); wire(edges(V, F), C_B, 3, "", "B")
    V, F = g["L"]; mesh(V, F, C_L, 0.12, "{L ≤ L(1/2)}  ⊇ B", "L"); wire(edges(V, F), C_L, 4, "", "L", dash="dash")
    # sliver and sorted curve
    V, F = g["sliver"]; mesh(V, F, C_SLIV, 0.22, "sorted sliver f₀ ≤ f₁ ≤ f₂ ≤ f₃", "sliv"); wire(edges(V, F), C_SLIV, 3, "", "sliv")
    S = g["sorted"]
    tr.append(go.Scatter3d(x=S[:, 0], y=S[:, 1], z=S[:, 2], mode="lines", line=dict(color="#333333", width=4, dash="dot"),
                           name="sort(f_p), the ordered version", legendgroup="sliv2", hoverinfo="skip"))
    # curve
    P, p = g["curve"], g["p"]
    tr.append(go.Scatter3d(x=P[:, 0], y=P[:, 1], z=P[:, 2], mode="lines",
                           line=dict(color=p, colorscale="Viridis", cmin=0, cmax=1, width=9),
                           name="binomial curve f_p (purple p=0 → yellow p=1)", legendgroup="curve",
                           text=[f"p = {v:.3f}" for v in p], hoverinfo="text"))
    pt = np.arange(0.1, 1.0, 0.1)
    Q = xyz(binomial(N, pt))
    tr.append(go.Scatter3d(x=Q[:, 0], y=Q[:, 1], z=Q[:, 2], mode="markers+text", text=[f"{v:.1f}" for v in pt],
                           textposition="top center", marker=dict(size=3, color="black"), showlegend=False,
                           legendgroup="curve", hoverinfo="skip"))
    # tie points
    for cusp in (False, True):
        ks = [k for k in TIES if (k in CUSPS) == cusp and k not in ((1, 2),)]
        Qt = xyz(binomial(N, [TIES[k] for k in ks]))
        lab = [("cusp " if cusp else "tie ") + (f"{k}, p = {TIES[k]:.3f}" if k != (0, 3) else "(0,3) and (1,2), p = 1/2")
               for k in ks]
        tr.append(go.Scatter3d(x=Qt[:, 0], y=Qt[:, 1], z=Qt[:, 2], mode="markers", hovertext=lab, hoverinfo="text",
                               marker=dict(size=8 if cusp else 6, symbol="circle-open" if cusp else "diamond-open",
                                           color="black", line=dict(width=3)),
                               name="cusps" if cusp else "tie points", legendgroup="ties"))
    # hidden layers
    for (i, j) in itertools.combinations(range(4), 2):
        W = wall(i, j)
        tr.append(go.Mesh3d(x=W[:, 0], y=W[:, 1], z=W[:, 2], i=[0], j=[1], k=[2], color="#777777", opacity=0.12,
                            name="walls f_i = f_j", legendgroup="walls", showlegend=(i, j) == (0, 1),
                            visible="legendonly", hoverinfo="name"))
    cc = ["#9b8fc7", "#7b6fb0", "#2a9d8f", "#e9a03b", "#d17a22", "#3a7dc9", "#6fa8dc", "#9fc5e8"]
    for (c, poly), col in zip(g["certs"].items(), cc):
        if poly is None: continue
        T = [(0, m, m + 1) for m in range(1, len(poly) - 1)]
        tr.append(go.Mesh3d(x=poly[:, 0], y=poly[:, 1], z=poly[:, 2], i=[t[0] for t in T], j=[t[1] for t in T],
                            k=[t[2] for t in T], color=col, opacity=0.3, name=f"certificate plane W_{c:g} = 2",
                            visible="legendonly", hoverinfo="name"))
    for t, (V, F) in g["levels"].items():
        wire(edges(V, F), "#888888", 2, f"level set E = {t:g}", f"lev{t}", show=False, legend=True)
    first = True
    for V, F in g["never"]:
        mesh(V, F, "#9aa0a6", 0.25, "16 non-unimodal chambers (never visited)", "never", show=False, legend=first)
        first = False

    fig = go.Figure(tr)
    ax = dict(visible=False)
    fig.update_layout(title="n = 3: the simplex of mass vectors (f₀, f₁, f₂, f₃) — drag to rotate, click the legend to show/hide",
                      scene=dict(xaxis=ax, yaxis=ax, zaxis=ax, aspectmode="data",
                                 camera=dict(eye=dict(x=0.25, y=-1.9, z=0.55))),
                      legend=dict(x=0.0, y=0.95, font=dict(size=13)), margin=dict(l=0, r=0, t=50, b=0),
                      paper_bgcolor="white")
    fig.write_html(path, include_plotlyjs=True)

# ------------------------------------------------------------------------------------------- matplotlib
def png(g, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
    plt.rcParams.update({"font.size": 16})
    views = [("from the front (curve's plane)", 8, -90), ("oblique", 24, -58), ("from above (down the symmetry axis)", 89, -90)]
    fig = plt.figure(figsize=(33, 12), dpi=100)
    fig.subplots_adjust(left=0, right=1, top=0.92, bottom=0.08, wspace=0)
    cmap = plt.get_cmap("viridis")
    for k, (title, el, az) in enumerate(views):
        ax = fig.add_subplot(1, 3, k + 1, projection="3d")
        ax.add_collection3d(Line3DCollection([TET[[i, j]] for i, j in itertools.combinations(range(4), 2)],
                                             colors="black", linewidths=1.6))
        V, F = g["L"]
        ax.add_collection3d(Poly3DCollection([V[f] for f in F], facecolors=C_L, alpha=0.07, edgecolors="none"))
        ax.add_collection3d(Line3DCollection(edges(V, F), colors=C_L, linewidths=1.4, linestyles="--"))
        V, F = g["B"]
        ax.add_collection3d(Poly3DCollection([V[f] for f in F], facecolors=C_B, alpha=0.22, edgecolors=C_B, linewidths=0.6))
        V, F = g["sliver"]
        ax.add_collection3d(Poly3DCollection([V[f] for f in F], facecolors=C_SLIV, alpha=0.18, edgecolors=C_SLIV, linewidths=1.0))
        S = g["sorted"]; ax.plot(S[:, 0], S[:, 1], S[:, 2], color="#333333", lw=1.6, ls=":")
        P, p = g["curve"], g["p"]
        lc = Line3DCollection(np.stack([P[:-1], P[1:]], 1), cmap=cmap, linewidths=4.5)
        lc.set_array(0.5*(p[:-1] + p[1:])); lc.set_clim(0, 1); ax.add_collection3d(lc)
        for kk, v in TIES.items():
            if kk == (1, 2): continue
            q = xyz(binomial(N, [v]))[0]
            ax.scatter(*q, s=110 if kk in CUSPS else 70, marker="o" if kk in CUSPS else "D",
                       facecolors="white", edgecolors="black", linewidths=2, depthshade=False, zorder=10)
        for i in range(4):
            ax.text(*(TET[i]*1.13), VNAME[i].split("  ")[0] + ("\np = 0" if i == 0 else "\np = 1" if i == 3 else ""),
                    ha="center", fontsize=15)
        ax.view_init(elev=el, azim=az)
        ax.set_box_aspect((1, 1, 0.75), zoom=1.22); ax.set_axis_off()
        lim = 1.05; ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-0.75, 0.75)
        ax.set_title(title, fontsize=19)
    from matplotlib.lines import Line2D
    hand = [Line2D([], [], color=cmap(0.5), lw=4.5, label="binomial curve f$_p$ (purple p = 0 → yellow p = 1)"),
            Line2D([], [], color=C_B, lw=8, alpha=0.4, label="B = {E ≤ E(1/2) = 2}"),
            Line2D([], [], color=C_L, lw=1.6, ls="--", label="{L ≤ L(1/2)} ⊇ B (edges)"),
            Line2D([], [], color=C_SLIV, lw=8, alpha=0.5, label="sorted sliver f₀ ≤ f₁ ≤ f₂ ≤ f₃"),
            Line2D([], [], color="#333333", lw=1.6, ls=":", label="sort(f$_p$)"),
            Line2D([], [], ls="", marker="D", ms=10, mfc="white", mec="black", mew=2, label="tie points"),
            Line2D([], [], ls="", marker="o", ms=12, mfc="white", mec="black", mew=2, label="cusps (1,3) and (0,2)")]
    fig.legend(handles=hand, loc="lower center", ncol=4, fontsize=16, frameon=False)
    fig.suptitle("n = 3: the simplex of mass vectors (f₀, f₁, f₂, f₃).  Rotate it in tetra_n3.html.",
                 x=0.02, ha="left", fontsize=22)
    fig.savefig(path, dpi=100, bbox_inches="tight")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    g = geometry()
    for name, fn in (("tetra_n3.html", html), ("tetra_n3.png", png)):
        path = os.path.join(a.out, name); fn(g, path); print(path)

if __name__ == "__main__":
    main()
