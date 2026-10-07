"""plotting/lower_bound/phantom_deficit.py -- what the phantom deficit is, and why the proof can afford it.

    .venv/bin/python plotting/lower_bound/phantom_deficit.py [--n 8] [--ns 8,20,50] [--out DIR]

A certificate with centre c (c in Z/2 + 1/4) ranks the positions 0..n by their distance from c.
  plain ranks  sigma_c(k): nearest real position gets n, the next n-1, ..., the farthest 0.
  formula      a_c(k) = n + 1/2 - 2|k - c|: ranks by distance on the WHOLE integer line, so integers
               outside 0..n ("phantoms") take rank slots too.  As a function of k it is a tent of
               slope 2 peaked at c.
While both sides of c still have real positions, the two agree: each step outward on one side is
matched by one on the other, and the ranks fall 2 per step.  Once one side runs out, every other slot
on that side goes to a phantom, the plain ranks on the far side fall only 1 per step, and the formula
falls behind:  sigma_c(k) - a_c(k) = #phantoms closer to c than k = max(0, t - k),  t = 2c - n - 1/2
(c > n/2; mirror image below).  Weighted by the masses this is the phantom deficit
    Delta_c(p) = sum_k (sigma_c(k) - a_c(k)) f_p(k) = E_p[(t - K)^+]  >= 0,
so W_c = <a_c, f_p> = <sigma_c, f_p> - Delta_c <= <sigma_c, f_p> <= E (the rearrangement freebie
survives: the formula only ever scores LOWER than the plain ranks).
Panels: (a) a central certificate -- no phantom is closer than any real position, no deficit;
(b) an off-centre one -- the tent and the true ranks part at the mirror point t; (c) the deficit
ramp weighted by f_p at a p where that certificate is active; (d) the deficit of L's active
certificate along p against L's margin L(p) - L(1/2), log scale: near 1/2, where the margin is
smallest, the deficit is ~2^-n.  Descriptive, double precision.
"""
import argparse, math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import binom

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _data import OUT

C_TENT, C_TRUE, C_PH, C_DEF = "#d1495b", "#111111", "#9aa0a6", "#e9a03b"

def sigma(c, n):
    o = np.argsort(np.abs(np.arange(n + 1) - c))
    r = np.empty(n + 1); r[o] = n - np.arange(n + 1)
    return r

def a_of(c, k, n):
    return n + 0.5 - 2.0*np.abs(np.asarray(k) - c)

def deficit_vec(c, n):
    return sigma(c, n) - a_of(c, np.arange(n + 1), n)

def certs(n):
    return np.arange(-1, 2*n + 1)/2.0 + 0.25

def active(n, p):
    """(L(p), active c, deficit of the active certificate) for an array of p."""
    k = np.arange(n + 1)
    cs = certs(n)
    A = np.array([a_of(c, k, n) for c in cs])
    S = np.array([sigma(c, n) for c in cs])
    f = binom.pmf(k[None, :], n, np.asarray(p)[:, None])
    W = f @ A.T
    j = np.argmax(W, axis=1)
    return W[np.arange(len(p)), j], cs[j], np.einsum("ij,ij->i", f, (S - A)[j])

def ladder(ax, n, c, title, ext=4):
    """Ranks against position on the integer line, phantoms included."""
    ks = np.arange(-ext, n + ext + 1)
    real = (ks >= 0) & (ks <= n)
    a = a_of(c, ks, n)
    xs = np.linspace(-ext - 0.4, n + ext + 0.4, 500)
    ax.plot(xs, a_of(c, xs, n), color=C_TENT, lw=1.6, alpha=0.6, zorder=1)
    ax.axvspan(-ext - 0.5, -0.5, color="#eef0f2", zorder=0); ax.axvspan(n + 0.5, n + ext + 0.5, color="#eef0f2", zorder=0)
    sg = sigma(c, n)
    d = deficit_vec(c, n)
    # deficit bars
    for k in range(n + 1):
        if d[k] > 0:
            ax.plot([k, k], [a[real][k], sg[k]], color=C_DEF, lw=6, alpha=0.75, solid_capstyle="butt", zorder=2)
    ax.plot(ks[~real], a[~real], "o", ms=13, mfc="white", mec=C_PH, mew=2, zorder=3)
    ax.plot(ks[real], a[real], "o", ms=13, color=C_TENT, zorder=4)
    ax.plot(np.arange(n + 1), sg, "s", ms=13, mfc="none", mec=C_TRUE, mew=2.4, zorder=5)
    # the slot numbers on the phantoms that matter (closer to c than some real position)
    far = sg.min()
    for kk, av in zip(ks[~real], a[~real]):
        if av > a_of(c, np.arange(n + 1), n).min() - 0.01:
            ax.annotate(f"{av:g}", (kk, av), xytext=(0, 13), textcoords="offset points", ha="center",
                        fontsize=14, color="#6b7280")
    ax.axvline(c, color=C_TENT, lw=1.2, ls=":")
    ax.text(c + 0.12, n + 1.05, f"c = {c:g}", ha="left", color=C_TENT, fontsize=16)
    ax.axhline(0, color="#bbbbbb", lw=1)
    ax.set_xticks(ks); ax.set_xlim(-ext - 0.5, n + ext + 0.5)
    ax.set_ylim(min(a.min(), -1) - 0.6, n + 1.6)
    ax.set_xlabel("position on the integer line (grey bands: phantoms, outside 0..n)")
    ax.set_title(title, loc="left", fontsize=19)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=8)
    ap.add_argument("--ns", default="8,20,50")
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    n = a.n
    c_mid, c_off = n/2 + 0.25, n - 1.25
    t_off = 2*c_off - n - 0.5

    plt.rcParams.update({"font.size": 16})
    fig = plt.figure(figsize=(32, 24), dpi=100)
    gs = fig.add_gridspec(2, 2, hspace=0.25, wspace=0.22)

    ax = fig.add_subplot(gs[0, 0])
    ladder(ax, n, c_mid, f"(a) a central certificate, c = n/2 + 1/4: plain ranks = formula, no deficit")
    ax.set_ylabel("rank / score given to the position")
    ax.annotate("every real position is closer to c\nthan any phantom: the tent IS the ranking",
                (round(c_mid) + 2, a_of(c_mid, round(c_mid) + 2, n)), xytext=(-4.2, n - 0.5), textcoords="data",
                fontsize=16, arrowprops=dict(arrowstyle="->", lw=1.4))

    bx = fig.add_subplot(gs[0, 1])
    ladder(bx, n, c_off, f"(b) an off-centre certificate, c = {c_off:g}: phantoms beyond n take slots")
    sg = sigma(c_off, n)
    bx.axvline(t_off, color=C_DEF, lw=1.6, ls="--")
    bx.text(t_off - 0.12, n + 1.05, f"t = 2c − n − ½ = {t_off:g}", ha="right", color=C_DEF, fontsize=16)
    ph = n + 1
    bx.annotate(f"phantom {ph} is closer to c than position {int(t_off) - 1}:\nit takes slot {a_of(c_off, ph, n):g}, "
                f"pushing {int(t_off) - 1} down to {a_of(c_off, int(t_off) - 1, n):g}",
                (ph, a_of(c_off, ph, n)), xytext=(n + 4.3, -7.5), textcoords="data", ha="right", fontsize=16,
                arrowprops=dict(arrowstyle="->", lw=1.4))
    bx.annotate("past t the formula keeps\nfalling 2 per step, the plain\nranks only 1: the gap (orange)\n"
                "grows by one per step,\nσ$_c$(k) − a$_c$(k) = max(0, t − k)",
                (1, (sg[1] + a_of(c_off, 1, n))/2), xytext=(-4.3, 1.0), textcoords="data", fontsize=16,
                arrowprops=dict(arrowstyle="->", lw=1.4))

    # (c) weighted by the masses
    cx = fig.add_subplot(gs[1, 0])
    k = np.arange(n + 1)
    ps = np.linspace(0.5, 0.999, 4000)
    _, cact, _ = active(n, ps)
    sel = ps[np.isclose(cact, c_off)]
    p_show = float(np.median(sel)) if len(sel) else c_off/n
    f = binom.pmf(k, n, p_show)
    d = deficit_vec(c_off, n)
    cx.bar(k - 0.2, f, width=0.38, color="#3a7dc9", alpha=0.75, label=f"masses f$_p$(k), p = {p_show:.3f} (where c = {c_off:g} is L's active centre)")
    cx2 = cx.twinx()
    cx2.bar(k + 0.2, d, width=0.38, color=C_DEF, alpha=0.75, label="deficit per position σ$_c$(k) − a$_c$(k)")
    cx.bar(k - 0.2, d*f, width=0.38, color="#7a1f2b", alpha=0.9, label="their product: each position's share of Δ")
    Dval = float(d @ f)
    cx.set_xticks(k); cx.set_xlabel("position k")
    cx.set_ylabel("mass"); cx2.set_ylabel("deficit (rank slots lost)")
    h1, l1 = cx.get_legend_handles_labels(); h2, l2 = cx2.get_legend_handles_labels()
    cx.legend(h1 + h2, l1 + l2, loc="upper left", fontsize=15, frameon=False)
    cx.set_title(f"(c) weighted by the masses: Δ = Σ (σ$_c$ − a$_c$) f$_p$ = E[(t − K)$^+$] = {Dval:.3f}\n"
                 f"the ramp sits in the far tail, opposite the bulk of the masses", loc="left", fontsize=19)
    cx.set_ylim(0, f.max()*1.45); cx2.set_ylim(0, d.max()*1.45 if d.max() > 0 else 1)

    # (d) along p, against the margin
    dx = fig.add_subplot(gs[1, 1])
    ns = [int(v) for v in a.ns.split(",")]
    cols = plt.get_cmap("viridis")(np.linspace(0.1, 0.8, len(ns)))
    for nn, col in zip(ns, cols):
        pp = 0.5 + np.linspace(0, 0.499, 6000)
        L, _, De = active(nn, pp)
        Lh = active(nn, np.array([0.5]))[0][0]
        x = pp - 0.5
        dx.plot(x, np.maximum(L - Lh, 1e-300), color=col, lw=2.4)
        dx.plot(x, np.maximum(De, 1e-300), color=col, lw=2.0, ls="--")
    dx.set_yscale("log"); dx.set_ylim(1e-30, 3)
    dx.set_xlim(0, 0.5)
    dx.set_xlabel("p − 1/2")
    dx.set_ylabel("size")
    hand = [Line2D([], [], color=col, lw=2.4, label=f"n = {nn}") for nn, col in zip(ns, cols)] + \
           [Line2D([], [], color="#555555", lw=2.4, label="L's margin L(p) − L(1/2)  (what (★) needs ≥ 0)"),
            Line2D([], [], color="#555555", lw=2.0, ls="--", label="phantom deficit of L's active certificate (what L gave up)")]
    dx.legend(handles=hand, loc="lower right", fontsize=15, frameon=False)
    dx.set_title("(d) Δ of L's active certificate against L's margin: exactly 0 at p = 1/2,\n"
                 "~2$^{-n}$ near it, where the margin is smallest; it grows only where the margin is large", loc="left", fontsize=19)

    hand = [Line2D([], [], color=C_TENT, lw=1.6, alpha=0.6, label="formula a$_c$(k) = n + ½ − 2|k − c|: a tent, slope 2"),
            Line2D([], [], ls="", marker="o", ms=12, color=C_TENT, label="formula score of a real position"),
            Line2D([], [], ls="", marker="o", ms=12, mfc="white", mec=C_PH, mew=2, label="phantom (outside 0..n): takes a slot too"),
            Line2D([], [], ls="", marker="s", ms=12, mfc="none", mec=C_TRUE, mew=2.2, label="plain rank σ$_c$(k) among the real positions"),
            Line2D([], [], color=C_DEF, lw=6, alpha=0.75, label="phantom deficit σ$_c$(k) − a$_c$(k) ≥ 0")]
    fig.legend(handles=hand, loc="upper center", ncol=3, fontsize=16, frameon=False, bbox_to_anchor=(0.5, 0.965))
    fig.suptitle(f"The phantom deficit (n = {n} in a–c): the formula ranks by distance on the whole integer line, so "
                 "positions beyond 0..n steal rank slots.\nIt only ever scores LOWER than the plain ranks, so L ≤ E is "
                 "untouched; it buys a closed form for each certificate, n + ½ − 2E|K − c|.",
                 x=0.01, y=1.0, ha="left", fontsize=21)
    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, f"phantom_deficit_n{n}.png")
    fig.savefig(path, dpi=100, bbox_inches="tight")
    print(path)

if __name__ == "__main__":
    main()
