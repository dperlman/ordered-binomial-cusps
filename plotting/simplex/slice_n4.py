"""plotting/simplex/slice_n4.py -- an exact 3D slice of the n = 4 simplex through the touch point p = 1/2.

    .venv/bin/python plotting/simplex/slice_n4.py [--n 4] [--near 0.1] [--out DIR]

The 4-simplex cannot be drawn whole, and projecting it would replace B by its shadow, which covers the
curve.  So it is cut by the 3-flat through x* = f_{1/2} spanned by the curve's tangent and the
symmetric directions (_slice.py).  That flat contains the curve's tangent and curvature at p = 1/2,
so near 1/2 the curve lies in it to second order.  Everything in the flat is exact: B, {L <= L(1/2)},
the walls, the simplex itself.  The curve is drawn as its orthogonal projection onto the flat, solid
and coloured by p while it is within --near (relative) of the flat, thin grey beyond that.
Axes: horizontal = along the curve's tangent at 1/2; vertical = away from the centre of the simplex
(the curve bends down, towards B); depth = the remaining symmetric direction.
The flat contains every symmetric vector, so ALL the walls f_k = f_{n-k} cut it in the same plane
(horizontal = 0): at p = 1/2 those pairs tie at once.  The sorted sliver only grazes the flat and is
not drawn.  Also runs for n = 5 (--n 5), where the flat is again 3D.
Writes simplex_figures/slice_nN.html (interactive) and slice_nN.png (three views + how far the
curve is from the flat).
"""
import argparse, csv, itertools, math, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _simplex import OUT, ROOT, binomial, certificates, perms, E_of
from _simplex3d import edges, triangles
from _slice import Slice, section

C_B, C_L = "#3a7dc9", "#d1495b"
SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")

def ties(n):
    """{(i, j): p*} for every pair, both sides of 1/2."""
    lnC = [math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1) for k in range(n + 1)]
    return {(i, j): 1/(1 + math.exp(-(lnC[i] - lnC[j])/(j - i))) for i, j in itertools.combinations(range(n + 1), 2)}

def cusps(n):
    with open(os.path.join(ROOT, "cusps", f"n{n:05d}.csv")) as fh:
        rows = [(int(r["i"]), int(r["j"]), float(r["pstar"])) for r in csv.DictReader(fh)]
    return {(i, j): p for i, j, p in rows} | {(n - j, n - i): 1 - p for i, j, p in rows}

def geometry(n, near):
    S = Slice(n)
    E0 = E_of(S.xs[None])[0]
    g = dict(S=S, E0=E0)
    g["frame"] = S.polytope([])
    g["B"] = S.polytope([(s, E0) for s in perms(n)])
    g["L"] = S.polytope([(a, E0) for a in certificates(n).values()])
    g["levels"] = {t: S.polytope([(s, t) for s in perms(n)]) for t in (E0 + 0.2, E0 + 0.4)}
    g["certs"] = {c: S.plane_in_slice(a, E0) for c, a in certificates(n).items()}
    walls = {}
    for i, j in itertools.combinations(range(n + 1), 2):
        a = np.zeros(n + 1); a[i], a[j] = 1, -1
        P = S.plane_in_slice(a, 0.0)
        if P is None: continue
        key = tuple(np.round(np.sort(P, axis=0).ravel(), 9))
        walls.setdefault(key, (P, []))[1].append((i, j))
    g["walls"] = list(walls.values())
    p = np.linspace(0, 1, 2001)
    f = binomial(n, p)
    g["p"], g["curve"], g["off"] = p, S.u(f), S.off(f)
    g["near"] = g["off"] <= near
    d = f - S.xs
    pl = S.A[:, [0, 2]]
    g["off2"] = np.linalg.norm(d - (d @ pl) @ pl.T, axis=1)/np.maximum(np.linalg.norm(d, axis=1), 1e-300)
    g["near2"] = g["off2"] <= near
    g["levels2"] = {t: S.polytope([(s_, t) for s_ in perms(n)]) for t in E0 + np.arange(0.1, 1.3, 0.1)}
    T, C = ties(n), cusps(n)
    g["ties"] = {k: (v, S.u(binomial(n, [v]))[0], S.off(binomial(n, [v]))[0]) for k, v in T.items()
                 if abs(v - 0.5) > 1e-12 and k not in C}             # cusps are drawn separately
    g["cusps"] = {k: (v, S.u(binomial(n, [v]))[0], S.off(binomial(n, [v]))[0]) for k, v in C.items()}
    g["axis_pairs"] = [k for k, v in T.items() if abs(v - 0.5) <= 1e-12]
    return g

def pair_str(k):
    return f"({k[0]},{k[1]})"

def wall_name(pairs):
    return ", ".join(f"f{i} = f{j}".translate(SUB) for i, j in pairs)

# ------------------------------------------------------------------------------------------- plotly
def html(g, n, path):
    import plotly.graph_objects as go
    tr = []
    def wire(segl, color, width, name, group, show=True, legend=False, dash="solid"):
        xs, ys, zs = [], [], []
        for s in segl:
            xs += [s[0][0], s[1][0], None]; ys += [s[0][1], s[1][1], None]; zs += [s[0][2], s[1][2], None]
        tr.append(go.Scatter3d(x=xs, y=ys, z=zs, mode="lines", line=dict(color=color, width=width, dash=dash),
                               name=name, legendgroup=group, showlegend=legend,
                               visible=True if show else "legendonly", hoverinfo="skip"))
    def mesh(V, F, color, opacity, name, group, show=True, legend=True):
        T = np.array(triangles(F))
        tr.append(go.Mesh3d(x=V[:, 0], y=V[:, 1], z=V[:, 2], i=T[:, 0], j=T[:, 1], k=T[:, 2], color=color,
                            opacity=opacity, flatshading=True, name=name, legendgroup=group, showlegend=legend,
                            visible=True if show else "legendonly", hoverinfo="name"))
    def poly(P, color, opacity, name, group, show, legend):
        T = [(0, m, m + 1) for m in range(1, len(P) - 1)]
        tr.append(go.Mesh3d(x=P[:, 0], y=P[:, 1], z=P[:, 2], i=[t[0] for t in T], j=[t[1] for t in T],
                            k=[t[2] for t in T], color=color, opacity=opacity, name=name, legendgroup=group,
                            showlegend=legend, visible=True if show else "legendonly", hoverinfo="name"))

    V, F = g["frame"]; wire(edges(V, F), "black", 4, f"the {n}-simplex, cut by the slice", "frame", legend=True)
    V, F = g["B"]; mesh(V, F, C_B, 0.45, "B = {E ≤ E(1/2)} in the slice", "B"); wire(edges(V, F), C_B, 2, "", "B")
    V, F = g["L"]; mesh(V, F, C_L, 0.12, "{L ≤ L(1/2)} in the slice  ⊇ B", "L"); wire(edges(V, F), C_L, 4, "", "L", dash="dash")
    # the curve
    P, p, off, nr = g["curve"], g["p"], g["off"], g["near"]
    hov = [f"p = {a:.3f}<br>off the slice: {100*b:.1f}%" for a, b in zip(p, off)]
    Pf = P.copy(); Pf[nr] = np.nan
    tr.append(go.Scatter3d(x=Pf[:, 0], y=Pf[:, 1], z=Pf[:, 2], mode="lines", line=dict(color="#9aa0a6", width=3),
                           name="curve, projected, where it is far from the slice", text=hov, hoverinfo="text"))
    Pn = P.copy(); Pn[~nr] = np.nan
    tr.append(go.Scatter3d(x=Pn[:, 0], y=Pn[:, 1], z=Pn[:, 2], mode="lines",
                           line=dict(color=p, colorscale="Viridis", cmin=0, cmax=1, width=9),
                           name="binomial curve, projected (purple p=0 → yellow p=1)", text=hov, hoverinfo="text"))
    pt = np.round(np.arange(0.05, 1.0, 0.05), 2)
    idx = [int(np.argmin(np.abs(p - v))) for v in pt]
    keep = [k for k in idx if nr[k]]
    tr.append(go.Scatter3d(x=P[keep, 0], y=P[keep, 1], z=P[keep, 2], mode="markers+text",
                           text=[f"{p[k]:.2f}" for k in keep], textposition="top center",
                           marker=dict(size=3, color="black"), showlegend=False, hoverinfo="skip"))
    # touch point, ties, cusps
    tr.append(go.Scatter3d(x=[0], y=[0], z=[0], mode="markers", marker=dict(size=7, color="black", symbol="square"),
                           name="p = 1/2: " + ", ".join(map(pair_str, g["axis_pairs"])) + " all tie; the touch point",
                           hovertext="x* = f_{1/2}", hoverinfo="text"))
    for name, d, sym, size in (("tie points", g["ties"], "diamond-open", 6), ("cusps", g["cusps"], "circle-open", 9)):
        ks = [k for k in d if d[k][2] <= 0.25]
        if not ks: continue
        Q = np.array([d[k][1] for k in ks])
        tr.append(go.Scatter3d(x=Q[:, 0], y=Q[:, 1], z=Q[:, 2], mode="markers", name=name,
                               marker=dict(size=size, symbol=sym, color="black", line=dict(width=3)),
                               hovertext=[f"{pair_str(k)}, p = {d[k][0]:.4f}, off the slice {100*d[k][2]:.1f}%" for k in ks],
                               hoverinfo="text"))
    # walls
    first = True
    for P_, pairs in g["walls"]:
        axis_wall = all(i + j == n for i, j in pairs)
        nm = wall_name(pairs) + ("  (all tie at p = 1/2)" if axis_wall else "")
        if axis_wall:
            poly(P_, "#555555", 0.10, "wall " + nm, "axiswall", True, True)
        else:
            poly(P_, "#777777", 0.10, "other walls f_i = f_j", "walls", False, first); first = False
    cc = ["#9b8fc7", "#7b6fb0", "#2a9d8f", "#e9a03b", "#d17a22", "#3a7dc9", "#6fa8dc", "#9fc5e8", "#b4a7d6", "#76a5af"]
    for (c, P_), col in zip(g["certs"].items(), cc):
        if P_ is not None:
            poly(P_, col, 0.3, f"certificate plane W_{c:g} = L(1/2)", f"c{c}", False, True)
    for t, (V, F) in g["levels"].items():
        wire(edges(V, F), "#888888", 2, f"level set E = {t:.4g}", f"lev{t}", show=False, legend=True)

    fig = go.Figure(tr)
    ax = dict(visible=False)
    fig.update_layout(title=(f"n = {n}: an exact 3D slice of the simplex through the touch point p = 1/2<br>"
                             "<sup>horizontal: along the curve; vertical: away from the centre; depth: the other "
                             "symmetric direction.  Drag to rotate; click the legend to show/hide.</sup>"),
                      title_font=dict(size=14),
                      scene=dict(xaxis=ax, yaxis=ax, zaxis=ax, aspectmode="data",
                                 camera=dict(eye=dict(x=0.35, y=-1.8, z=0.45))),
                      legend=dict(x=0.0, y=0.95, font=dict(size=12)), margin=dict(l=0, r=0, t=50, b=0))
    fig.write_html(path, include_plotlyjs=True)

# ------------------------------------------------------------------------------------------- matplotlib
def seg_at_depth(P):
    """Where a planar polygon P (k, 3) in the slice crosses depth u[1] = 0: a segment in (u0, u2)."""
    if P is None: return None
    s = P[:, 1]; pts = []
    for i in range(len(P)):
        j = (i + 1) % len(P)
        if s[i]*s[j] < 0: pts.append(P[i] + (P[j] - P[i])*s[i]/(s[i] - s[j]))
        elif abs(s[i]) < 1e-12: pts.append(P[i])
    if len(pts) < 2: return None
    pts = np.array(pts)[:, [0, 2]]
    d = ((pts[:, None] - pts[None])**2).sum(-1)
    i, j = np.unravel_index(np.argmax(d), d.shape)
    return pts[[i, j]]

def draw_section(ax, g, n, near, zoom=False):
    from matplotlib.patches import Polygon
    from matplotlib.collections import LineCollection
    import matplotlib.pyplot as plt
    cmap = plt.get_cmap("viridis")
    sec = lambda VF: section(*VF, axis=1)
    for t, VF in sorted(g["levels2"].items()):
        Q = sec(VF)
        if Q is not None: ax.add_patch(Polygon(Q, closed=True, fc="none", ec="#a3a8ae", lw=1.0, zorder=1))
    Q = sec(g["L"]); ax.add_patch(Polygon(Q, closed=True, fc=C_L, alpha=0.10, ec="none", zorder=2))
    ax.add_patch(Polygon(Q, closed=True, fc="none", ec=C_L, lw=2.6, ls=(0, (6, 3)), zorder=5))
    Q = sec(g["B"]); ax.add_patch(Polygon(Q, closed=True, fc=C_B, alpha=0.28, ec=C_B, lw=2.6, zorder=3))
    Q = sec(g["frame"]); ax.add_patch(Polygon(Q, closed=True, fc="none", ec="black", lw=2.2, zorder=6))
    for P_, pairs in g["walls"]:
        sg = seg_at_depth(P_)
        if sg is None: continue
        axis_wall = all(i + j == n for i, j in pairs)
        ax.plot(sg[:, 0], sg[:, 1], color="#555555", lw=1.6 if axis_wall else 1.0,
                ls=(0, (2, 3)), alpha=0.9 if axis_wall else 0.5, zorder=4)
    for c, P_ in g["certs"].items():
        sg = seg_at_depth(P_)
        if sg is not None:
            ax.plot(sg[:, 0], sg[:, 1], color="#2a9d8f", lw=1.0, alpha=0.6, zorder=4)
    P2, p, nr = g["curve"][:, [0, 2]], g["p"], g["near2"]
    far = np.where(nr[:, None], np.nan, P2)
    ax.plot(far[:, 0], far[:, 1], color="#9aa0a6", lw=1.6, zorder=7)
    m = nr[:-1] & nr[1:]
    lc = LineCollection(np.stack([P2[:-1], P2[1:]], 1)[m], cmap=cmap, linewidths=5, zorder=8, capstyle="round")
    lc.set_array(0.5*(p[:-1] + p[1:])[m]); lc.set_clim(0, 1); ax.add_collection(lc)
    ax.plot(0, 0, "s", ms=11, color="black", zorder=11)
    for d, mk, sz in ((g["ties"], "D", 10), (g["cusps"], "o", 14)):
        for kk, (pv, u, o) in d.items():
            if g["off2"][int(np.argmin(np.abs(p - pv)))] <= near:
                ax.plot(u[0], u[2], mk, ms=sz, mfc="white", mec="black", mew=2.2, zorder=10)
                if zoom or d is g["cusps"]:
                    ax.annotate(f"{pair_str(kk)}  p = {pv:.3f}", (u[0], u[2]), xytext=(0, 16 if d is g["cusps"] else -24),
                                textcoords="offset points", ha="center", fontsize=15, zorder=12)
    for t in np.round(np.arange(0.35, 0.66, 0.05), 2):
        k = int(np.argmin(np.abs(p - t)))
        if nr[k] and abs(t - 0.5) > 1e-9 and zoom:
            ax.plot(*P2[k], "o", ms=4, color="black", zorder=9)
    ax.set_aspect("equal")

def png(g, n, near, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle
    plt.rcParams.update({"font.size": 16})
    cmap = plt.get_cmap("viridis")
    fig = plt.figure(figsize=(34, 26), dpi=100)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], hspace=0.12, wspace=0.08)
    P, p = g["curve"], g["p"]

    # (a) the osculating plane, whole
    ax = fig.add_subplot(gs[0, 0]); draw_section(ax, g, n, near)
    Q = section(*g["frame"], axis=1); lo, hi = Q.min(0), Q.max(0); pad = 0.04*(hi - lo).max()
    ax.set_xlim(lo[0] - pad, hi[0] + pad); ax.set_ylim(lo[1] - pad, hi[1] + pad)
    ax.set_title("(a) the curve's own plane (depth 0): every region cut exactly", loc="left", fontsize=20)
    ax.set_xlabel("along the curve's tangent at p = 1/2"); ax.set_ylabel("away from the centre")
    # (b) zoom on the touch point
    sel = (p > 0.36) & (p < 0.64)
    bx0, bx1 = P[sel, 0].min(), P[sel, 0].max(); bz0, bz1 = P[sel, 2].min(), P[sel, 2].max()
    w = bx1 - bx0; zx0, zx1 = bx0 - 0.08*w, bx1 + 0.08*w
    hz = (zx1 - zx0)*0.62; zz0, zz1 = bz1 - 0.75*hz, bz1 + 0.25*hz
    ax.add_patch(Rectangle((zx0, zz0), zx1 - zx0, zz1 - zz0, fill=False, ec="#2b6cb0", lw=2.2, zorder=12))
    bx = fig.add_subplot(gs[0, 1]); draw_section(bx, g, n, near, zoom=True)
    bx.set_xlim(zx0, zx1); bx.set_ylim(zz0, zz1)
    bx.set_title("(b) close up on the touch point (blue box in a); dots every 0.05 in p", loc="left", fontsize=20)
    bx.set_xlabel("along the curve's tangent at p = 1/2")

    # (c) the 3D slice, oblique
    cx = fig.add_subplot(gs[1, 0], projection="3d")
    V0, F0 = g["frame"]; lo3, hi3 = V0.min(0), V0.max(0)
    cx.add_collection3d(Line3DCollection(edges(V0, F0), colors="black", linewidths=1.4))
    V, F = g["L"]; cx.add_collection3d(Line3DCollection(edges(V, F), colors=C_L, linewidths=1.2, linestyles="--"))
    V, F = g["B"]; cx.add_collection3d(Poly3DCollection([V[f] for f in F], facecolors=C_B, alpha=0.22, edgecolors=C_B, linewidths=0.5))
    nr = g["near"]; m = nr[:-1] & nr[1:]
    lc = Line3DCollection(np.stack([P[:-1], P[1:]], 1)[m], cmap=cmap, linewidths=4.5)
    lc.set_array(0.5*(p[:-1] + p[1:])[m]); lc.set_clim(0, 1); cx.add_collection3d(lc)
    cx.scatter(0, 0, 0, marker="s", s=60, color="black", depthshade=False)
    cx.view_init(elev=22, azim=-60)
    cx.set_xlim(lo3[0], hi3[0]); cx.set_ylim(lo3[1], hi3[1]); cx.set_zlim(lo3[2], hi3[2])
    cx.set_box_aspect(tuple(hi3 - lo3), zoom=1.1); cx.set_axis_off()
    cx.set_title("(c) the whole 3D slice (rotate it in the HTML); (a) is its cut at depth 0", loc="left", fontsize=20)

    # (d) how far the curve is from the slice and from the plane
    dx = fig.add_subplot(gs[1, 1])
    dx.plot(p, 100*g["off"], color="#333333", lw=2.2, label="from the 3D slice (c)")
    dx.plot(p, 100*g["off2"], color="#333333", lw=1.4, ls="--", label="from the plane (a), (b)")
    dx.axhline(100*near, color="#9aa0a6", ls=":", lw=1.4)
    for d, mk, sz in ((g["ties"], "D", 9), (g["cusps"], "o", 13)):
        for kk, (pv, u, o) in d.items():
            dx.plot(pv, 100*o, mk, ms=sz, mfc="white", mec="black", mew=2)
    dx.set_xlim(0, 1); dx.set_ylim(0, 105*max(g["off"].max(), g["off2"].max()))
    dx.set_xlabel("p"); dx.set_ylabel("distance of f$_p$ off the cut, % of |f$_p$ − x*|")
    dx.legend(loc="upper center", fontsize=16, frameon=False)
    dx.set_title(f"(d) how faithful the picture is: solid curve while under {100*near:.0f}%, grey beyond", loc="left", fontsize=20)

    hand = [Line2D([], [], color=cmap(0.5), lw=5, label="binomial curve, projected (purple p = 0 → yellow p = 1)"),
            Line2D([], [], color=C_B, lw=8, alpha=0.4, label="B = {E ≤ E(1/2)}"),
            Line2D([], [], color=C_L, lw=2.4, ls="--", label="{L ≤ L(1/2)}"),
            Line2D([], [], color="#a3a8ae", lw=1.2, label="level sets of E"),
            Line2D([], [], color="#2a9d8f", lw=1.2, alpha=0.7, label="certificate lines W$_c$ = L(1/2)"),
            Line2D([], [], color="#555555", lw=1.4, ls=(0, (2, 3)), label="walls f$_i$ = f$_j$ (darker: all f$_k$ = f$_{n-k}$ at once)"),
            Line2D([], [], color="black", lw=2.2, label="the simplex, cut"),
            Line2D([], [], ls="", marker="s", ms=10, color="black", label="x* = f$_{1/2}$, the touch point"),
            Line2D([], [], ls="", marker="D", ms=10, mfc="white", mec="black", mew=2, label="tie points"),
            Line2D([], [], ls="", marker="o", ms=13, mfc="white", mec="black", mew=2, label="cusps")]
    fig.legend(handles=hand, loc="upper center", ncol=5, fontsize=16, frameon=False, bbox_to_anchor=(0.5, 0.955))
    fig.suptitle(f"n = {n}: the simplex cut through the touch point p = 1/2.  The cut contains the curve's tangent "
                 "and curvature there, so near 1/2 the curve lies in it to second order.",
                 x=0.01, y=0.99, ha="left", fontsize=22)
    fig.savefig(path, dpi=100, bbox_inches="tight")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=4)
    ap.add_argument("--near", type=float, default=0.10, help="draw the curve solid while off(p) <= this")
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    g = geometry(a.n, a.near)
    for name, fn in ((f"slice_n{a.n}.html", lambda g_, pth: html(g_, a.n, pth)),
                     (f"slice_n{a.n}.png", lambda g_, pth: png(g_, a.n, a.near, pth))):
        path = os.path.join(a.out, name); fn(g, path); print(path)

if __name__ == "__main__":
    main()
