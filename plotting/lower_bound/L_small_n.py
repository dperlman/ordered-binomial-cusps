"""plotting/lower_bound/L_small_n.py -- L(p) and the certificates it is built from, under E(n,p), at one small n.

    .venv/bin/python plotting/lower_bound/L_small_n.py [--n 12] [--zoom-xmax X] [--out DIR]

What is drawn (RESEARCH_LOG.md section 4 and the 2026-09-24 claude.ai entry, items 1-3):
  certificates         W_c(p) = n + 1/2 - 2 E_p|K - c|,  c in Z/2 + 1/4, c = -1/4 .. n+1/4.
                       Each one is a lower bound for E at every p (distance-from-c ranking,
                       rearrangement inequality); faint, coloured by c.
  L                    L(p) = max_c W_c(p) = n + 1/2 - 2 D(p), the certificate envelope;
                       L <= E everywhere, L(1/2) = E(1/2).
  switch points        p_m (star_check.switch_points), where the maximising c moves; L's kinks.
  E(n,p)               binom_core.E_at on the grid; cusps (local minima) from cusps/nNNNNN.csv,
                       every other tie point as a small tick.
Panels: left column the whole of [0,1], right column a close-up above 1/2 with a top axis in
x = n(p - 1/2) (the variable of L_scaled.py and L_arches.py).  Top row the curves, bottom row the gap E - L.
The grid is uniform plus every tie point and switch point, so every kink is drawn exactly.
Descriptive, double precision, not certified.
"""
import argparse, csv, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import binom

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _data import ROOT, OUT, core, switch_points

def tie_points(n):
    """p* of every pair 0 <= i < j <= n (both sides of 1/2; i+j = n gives exactly 1/2)."""
    l = core.lnC_arr(n)
    i, j = np.triu_indices(n + 1, 1)
    return 1.0/(1.0 + np.exp(-(l[i] - l[j])/(j - i)))

def cusps(n):
    """(p*, E) of the cusps above 1/2, from the certified per-n table."""
    path = os.path.join(ROOT, "cusps", f"n{n:05d}.csv")
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    return np.array([float(r["pstar"]) for r in rows]), np.array([float(r["E"]) for r in rows])

def family(n, p):
    """W_c(p) for every c, shape (len(c), len(p))."""
    c = np.arange(-1, 2*n + 1)/2.0 + 0.25
    k = np.arange(n + 1)
    f = binom.pmf(k[None, :], n, p[:, None])                   # (P, n+1)
    mad = np.abs(k[None, :] - c[:, None]) @ f.T                 # (C, P)
    return c, n + 0.5 - 2.0*mad

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=12)
    ap.add_argument("--zoom-xmax", type=float, default=None,
                    help="close-up (sampled at --grid points of its own) up to x = n(p-1/2) (default: just past the last cusp)")
    ap.add_argument("--grid", type=int, default=6000)
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    n = a.n

    # ---- data
    sw = switch_points(n)[3]
    sw_all = np.concatenate([sw, 1.0 - sw])
    ties = tie_points(n)
    cp, cE = cusps(n)
    zx = a.zoom_xmax if a.zoom_xmax else n*(cp.max() - 0.5) + 0.75
    zhi = min(1.0, 0.5 + zx/n)
    p = np.unique(np.concatenate([np.linspace(0, 1, a.grid + 1)[1:-1], np.linspace(0.5, zhi, a.grid),
                                  ties, sw_all, [0.5]]))
    p = p[(p > 0) & (p < 1)]
    E = np.array([core.E_at(n, v) for v in p])
    c, W = family(n, p)
    L = W.max(axis=0)
    Eh = core.E_half(n)
    # endpoints: E(0) = E(1) = n; D(0) = D(1) = 1/4, so L = n there too
    p = np.concatenate([[0.0], p, [1.0]])
    E = np.concatenate([[n], E, [n]])
    W = np.concatenate([family(n, np.array([1e-300]))[1], W, family(n, np.array([1 - 1e-16]))[1]], axis=1)
    L = W.max(axis=0)
    cp_all = np.concatenate([cp, 1 - cp]); cE_all = np.concatenate([cE, cE])
    Et = np.interp(ties, p, E)
    Lsw = np.interp(sw_all, p, L)

    # ---- figure
    plt.rcParams.update({"font.size": 20})
    fig, ax = plt.subplots(2, 2, figsize=(32, 20), dpi=100,
                           gridspec_kw=dict(height_ratios=[2.2, 1], hspace=0.28, wspace=0.14))
    cmap = plt.get_cmap("viridis")
    col = lambda ci: cmap(ci/(len(c) - 1))
    CE, CL = "#111111", "#d1495b"

    def curves(axh, shift, lo, hi):
        m = (p >= lo) & (p <= hi)
        for ci in range(len(c)):
            axh.plot(p[m], W[ci, m] - shift, color=col(ci), lw=1.6, alpha=0.45, zorder=1)
        axh.plot(p[m], E[m] - shift, color=CE, lw=3.2, zorder=4)
        axh.plot(p[m], L[m] - shift, color=CL, lw=3.2, zorder=5)
        for s in sw_all[(sw_all >= lo) & (sw_all <= hi)]:
            axh.axvline(s, color=CL, lw=1.0, ls=":", alpha=0.7, zorder=0)
        mt = (ties >= lo) & (ties <= hi)
        axh.plot(ties[mt], Et[mt] - shift, "|", color=CE, ms=14, mew=1.6, zorder=6)
        mc = (cp_all >= lo) & (cp_all <= hi)
        axh.plot(cp_all[mc], cE_all[mc] - shift, "o", ms=13, mfc="white", mec=CE, mew=2.5, zorder=7)
        ms = (sw_all >= lo) & (sw_all <= hi)
        axh.plot(sw_all[ms], Lsw[ms] - shift, "D", ms=9, color=CL, zorder=7)
        axh.axhline(0 if shift else Eh, color="#888888", lw=1.2, ls="--", zorder=0)
        axh.set_xlim(lo, hi)

    # A: everything
    A = ax[0, 0]
    curves(A, 0.0, 0.0, 1.0)
    A.set_ylim(Eh - 0.45*(n - Eh), n + 0.12*(n - Eh))
    A.set_ylabel("value")
    A.set_title(f"n = {n}: E(n,p), L(p) = max$_c$ W$_c$(p), and the certificates W$_c$", loc="left")

    # B: close-up above 1/2, relative to E(1/2)
    B = ax[0, 1]
    curves(B, Eh, 0.5, zhi)
    m = (p >= 0.5) & (p <= zhi)
    top = max((E[m] - Eh).max(), 1e-12)
    B.set_ylim(-0.6*top, 1.08*top)
    B.set_ylabel("value $-$ E(n,1/2)")
    B.set_title("close-up above 1/2 (dotted: switch points of L)", loc="left", pad=70)
    sec = B.secondary_xaxis("top", functions=(lambda v: n*(v - 0.5), lambda x: 0.5 + x/n))
    sec.set_xlabel("x = n(p $-$ 1/2)")
    for k, s in enumerate(sw[sw <= zhi], 1):
        B.annotate(f"m={k}", (s, 1.08*top), xytext=(4, -26), textcoords="offset points",
                   color=CL, fontsize=17)

    # C, D: the gap E - L
    for axh, lo, hi in ((ax[1, 0], 0.0, 1.0), (ax[1, 1], 0.5, zhi)):
        mm = (p >= lo) & (p <= hi)
        axh.fill_between(p[mm], 0, E[mm] - L[mm], color="#9aa5b1", alpha=0.5, lw=0)
        axh.plot(p[mm], E[mm] - L[mm], color="#3d434b", lw=2.2)
        for s in sw_all[(sw_all >= lo) & (sw_all <= hi)]:
            axh.axvline(s, color=CL, lw=1.0, ls=":", alpha=0.7)
        mc = (cp_all >= lo) & (cp_all <= hi)
        axh.plot(cp_all[mc], np.interp(cp_all[mc], p, E - L), "o", ms=11, mfc="white", mec=CE, mew=2.2)
        axh.set_xlim(lo, hi); axh.set_ylim(bottom=0)
        axh.set_xlabel("p")
        axh.axvline(0.5, color="#888888", lw=1.2, ls="--")
    ax[1, 0].set_ylabel("E $-$ L  (slack of L)")
    ax[1, 1].set_ylim(top=1.1*(E[m] - L[m]).max())
    ax[0, 1].sharex(ax[1, 1])
    ax[0, 0].sharex(ax[1, 0])

    hand = [Line2D([], [], color=CE, lw=3.2, label="E(n,p)"),
            Line2D([], [], color=CL, lw=3.2, label="L(p) = n + 1/2 $-$ 2 D(p), the certificate envelope"),
            Line2D([], [], color=col(len(c)//2), lw=1.6, alpha=0.6,
                   label=f"certificates W$_c$, c = $-$1/4 .. {n}+1/4 (colour by c)"),
            Line2D([], [], ls="", marker="o", ms=13, mfc="white", mec=CE, mew=2.5, label="cusps of E"),
            Line2D([], [], ls="", marker="|", ms=14, mew=1.6, color=CE, label="other tie points of E"),
            Line2D([], [], ls="", marker="D", ms=9, color=CL, label="switch points of L (kinks)"),
            Line2D([], [], color="#888888", lw=1.2, ls="--", label="E(n,1/2) = L(1/2)")]
    A.legend(handles=hand, loc="upper center", fontsize=18, framealpha=0.92)

    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, f"L_small_n_n{n:03d}.png")
    fig.savefig(path, dpi=100, bbox_inches="tight")
    print(path)

if __name__ == "__main__":
    main()
