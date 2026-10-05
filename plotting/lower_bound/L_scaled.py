"""plotting/lower_bound/L_scaled.py -- E and L near 1/2 in the trough scaling n^1.5.

    .venv/bin/python plotting/lower_bound/L_scaled.py [--ns 30,31,100,101,1000,1001]
                     [--xmax 3] [--troughs 1,2] [--umax 2] [--grid 3000] [--workers 8] [--out DIR]

L = n + 1/2 - 2 D(p) is the certificate envelope of the (★) reduction; data from _data.compute.
Coordinates (RESEARCH_LOG 2026-09-24, claude.ai entry item 4 and the Claude Code entry, result 4):
  x = n (p - 1/2)
  y_E = (E(n,p) - E(n,1/2)) n^1.5          the true margin
  y_L = (L(p) - L(1/2)) n^1.5 = 2 (D(1/2) - D(p)) n^1.5     L's margin;  y_L <= y_E
They do NOT collapse as curves: between switch points both form arches whose height in these
units grows ~n (L_arches.py shows them in their own scaling).  Only the TROUGHS (x ~ m/2, where the
cusps and switch points sit) settle onto the smooth term 2 sqrt(2/pi) x^2 (dotted).  At the m-th
switch point y_L is the scaled margin M(n,m) of star_check.py, whose limit is
  (m^2 - 1/4)/sqrt(2 pi)  [m odd, n even],  (m^2 + 1/4)/sqrt(2 pi)  [m odd, n odd],  m^2/sqrt(2 pi)  [m even]
-- drawn as black crosses.  Top row: x in [0, xmax], even n left, odd n right.  Bottom row: the
troughs listed in --troughs close up, in u = n (x - x_m) = n^2 (p - p_m) (the slopes at a trough
are ~n in x), y unchanged.  Output: lower_bound_proof_plots/.
"""
import argparse, math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _data import compute, trough_limit as limit, K2, OUT

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ns", default="30,31,100,101,1000,1001")
    ap.add_argument("--xmax", type=float, default=3.0)
    ap.add_argument("--grid", type=int, default=3000)
    ap.add_argument("--troughs", default="1,2", help="switch points m shown close up (bottom row)")
    ap.add_argument("--umax", type=float, default=2.0, help="half-width of a trough close-up in u")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    ns = sorted(int(v) for v in a.ns.split(","))

    S = {}
    for n, d in compute(ns, a.xmax, a.grid, local=a.umax, workers=a.workers).items():
        s = n**1.5
        S[n] = dict(x=d["x"], yE=d["dE"]*s, yL=d["dL"]*s, cx=d["cx"], cy=d["cdE"]*s,
                    m=d["m"], xm=d["xm"], ym=d["dLm"]*s)
        print(f"    n={n}: min y_L at a switch point = {S[n]['ym'].min():.6f} "
              f"(m={d['m'][np.argmin(S[n]['ym'])]})")

    # ---- figure: top row the whole range, bottom row each trough close up in u = n(x - x_m)
    troughs = [int(v) for v in a.troughs.split(",") if v]
    pars = (("even", [n for n in ns if n % 2 == 0]), ("odd", [n for n in ns if n % 2 == 1]))
    bottom = [f"{par}{m}" for m in troughs for par, _ in pars]
    w = len(bottom)
    mosaic = [["even"]*(w//2) + ["odd"]*(w - w//2), bottom]
    plt.rcParams.update({"font.size": 20})
    fig, ax = plt.subplot_mosaic(mosaic, figsize=(32, 22), dpi=100,
                                 gridspec_kw=dict(height_ratios=[1.5, 1], hspace=0.2, wspace=0.18))
    xs = np.linspace(0, a.xmax, 600)
    ramp = lambda name, k: [plt.get_cmap(name)(0.45 + 0.5*t/max(k - 1, 1)) for t in range(k)]

    for par, group in pars:
        top = ax[par]
        greys, reds = ramp("Greys", len(group)), ramp("Reds", len(group))
        n0 = group[0] if group else (0 if par == "even" else 1)
        mm = np.arange(1, int(2*a.xmax) + 1)
        top.plot(xs, K2*xs**2, ":", color="#2b6cb0", lw=2.5, zorder=1)
        top.plot(mm/2, limit(mm, n0), "x", color="black", ms=16, mew=3, zorder=8)
        for t, n in enumerate(group):
            d = S[n]
            top.plot(d["x"], d["yE"], color=greys[t], lw=2.2, zorder=3)
            top.plot(d["x"], d["yL"], color=reds[t], lw=2.2, zorder=4)
            top.plot(d["cx"], d["cy"], "o", ms=10, mfc="white", mec=greys[t], mew=2.2, zorder=6)
            top.plot(d["xm"], d["ym"], "D", ms=8, color=reds[t], zorder=7)
        top.set_xlim(0, a.xmax); top.set_ylim(-0.03*K2*a.xmax**2, 1.08*K2*a.xmax**2)
        top.axhline(0, color="#888888", lw=1)
        for h in np.arange(0.5, a.xmax + 0.01, 0.5):
            top.axvline(h, color="#cccccc", lw=1, zorder=0)
        top.set_xlabel("x = n(p $-$ 1/2)")
        top.set_title(f"{par} n: " + ", ".join(map(str, group)), loc="left")
        hand = [Line2D([], [], color=greys[t], lw=3, label=f"E,  n = {n}") for t, n in enumerate(group)] + \
               [Line2D([], [], color=reds[t], lw=3, label=f"L,  n = {n}") for t, n in enumerate(group)] + \
               [Line2D([], [], ls=":", color="#2b6cb0", lw=2.5, label="2$\\sqrt{2/\\pi}$ x$^2$"),
                Line2D([], [], ls="", marker="o", ms=10, mfc="white", mec="#555555", mew=2.2, label="cusps of E"),
                Line2D([], [], ls="", marker="D", ms=8, color="#c0392b", label="switch points of L"),
                Line2D([], [], ls="", marker="x", ms=16, mew=3, color="black",
                       label=f"limit of L at switch m: (m$^2${'$-$' if par == 'even' else '+'}[m odd]/4)/$\\sqrt{{2\\pi}}$")]
        top.legend(handles=hand, loc="upper left", fontsize=17, framealpha=0.92, ncol=2)

        for m in troughs:
            b = ax[f"{par}{m}"]
            lo, hi = np.inf, -np.inf
            for t, n in enumerate(group):
                d = S[n]
                if m > len(d["xm"]): continue
                xm = d["xm"][m - 1]
                u = n*(d["x"] - xm)
                k = np.abs(u) <= a.umax
                b.plot(u[k], d["yE"][k], color=greys[t], lw=2.4, zorder=3)
                b.plot(u[k], d["yL"][k], color=reds[t], lw=2.4, zorder=4)
                cu = n*(d["cx"] - xm); kc = np.abs(cu) <= a.umax
                b.plot(cu[kc], d["cy"][kc], "o", ms=11, mfc="white", mec=greys[t], mew=2.2, zorder=6)
                b.plot([0], [d["ym"][m - 1]], "D", ms=9, color=reds[t], zorder=7)
                lo = min(lo, d["yL"][k].min()); hi = max(hi, d["yE"][k].max())
            b.axhline(limit(m, n0), color="black", lw=1.5, ls="--", zorder=1)
            b.axhline(K2*(m/2)**2, color="#2b6cb0", lw=2.5, ls=":", zorder=1)
            b.axvline(0, color="#cccccc", lw=1, zorder=0)
            pad = 0.08*(hi - lo)
            b.set_ylim(min(lo, limit(m, n0)) - pad, hi + pad)
            b.set_xlim(-a.umax, a.umax)
            b.set_xlabel("u = n (x $-$ x$_m$)")
            b.set_title(f"trough m = {m}, {par} n  (dashed: limit of L)", loc="left", fontsize=19)
    ax["even"].set_ylabel("(value $-$ value(1/2)) n$^{1.5}$")
    ax[bottom[0]].set_ylabel("same, close up")

    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, f"L_scaled_x{a.xmax:g}_n{'-'.join(map(str, ns))}.png")
    fig.savefig(path, dpi=100, bbox_inches="tight")
    print(path)

if __name__ == "__main__":
    main()
