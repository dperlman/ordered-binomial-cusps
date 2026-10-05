"""plotting/lower_bound/L_arches.py -- the arches of E and L near 1/2, in their own scaling sqrt(n).

    .venv/bin/python plotting/lower_bound/L_arches.py [--ns 30,31,100,101,1000,1001,3000,3001]
                     [--xmax 3] [--grid 3000] [--workers 8] [--out DIR]

L = n + 1/2 - 2 D(p) is the certificate envelope of the (★) reduction; data from _data.compute.
x = n (p - 1/2).  Rows:
  (a) sqrt(n) (L - L(1/2)) against x.  Converges at fixed x to the periodic chain of parabolic arches
          A(x) = (1/4 - 4 d(x)^2)/sqrt(2 pi),   d(x) = distance from x to the nearest point of Z/2 + 1/4,
      height 1/(4 sqrt(2 pi)) = 0.0997 at x in Z/2 + 1/4, zero at x in Z/2 (the switch points).
      Heuristic: with the active centre c fixed, the local CLT gives E_p|K-c| ~ sigma sqrt(2/pi)
      + 2 (c-np)^2/sqrt(2 pi n) + (a lattice term that is the same for every quarter-point c), and
      L(p) - L(1/2) = 2 (D(1/2) - D(p)).  E coincides with L on most of each arch (E - L = 0
      exactly where the distance ranking is the true ranking), so its curve is not drawn here.
  (b) the next order for L:  n^1.5 (L - L(1/2)) - n A(x) - 2 sqrt(2/pi) x^2.
  (c) the same for E.  E - L is O(n^-1.5) near the troughs, so it shows up only at this order.
Even n in blues, odd n in oranges; darker = larger n.  Output: lower_bound_proof_plots/.
"""
import argparse, math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _data import compute, arch, K2, OUT

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ns", default="30,31,100,101,1000,1001,3000,3001")
    ap.add_argument("--xmax", type=float, default=3.0)
    ap.add_argument("--grid", type=int, default=3000)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    ns = sorted(int(v) for v in a.ns.split(","))
    D = compute(ns, a.xmax, a.grid, local=3.0, workers=a.workers)

    evens = [n for n in ns if n % 2 == 0]; odds = [n for n in ns if n % 2 == 1]
    ramp = lambda name, k: [plt.get_cmap(name)(0.4 + 0.55*t/max(k - 1, 1)) for t in range(k)]
    colour = {**dict(zip(evens, ramp("Blues", len(evens)))), **dict(zip(odds, ramp("Oranges", len(odds))))}

    plt.rcParams.update({"font.size": 20})
    fig, ax = plt.subplots(3, 1, figsize=(30, 26), dpi=100, sharex=True,
                           gridspec_kw=dict(height_ratios=[1.4, 1, 1], hspace=0.12))
    xs = np.linspace(0, a.xmax, 4001)
    ax[0].plot(xs, arch(xs), color="black", lw=2.2, ls="--", zorder=10)
    for n in ns:
        d = D[n]; x = d["x"]
        ax[0].plot(x, math.sqrt(n)*d["dL"], color=colour[n], lw=2.4)
        base = n*arch(x) + K2*x**2
        ax[1].plot(x, n**1.5*d["dL"] - base, color=colour[n], lw=2.0)
        ax[2].plot(x, n**1.5*d["dE"] - base, color=colour[n], lw=2.0)
        cb = n*arch(d["cx"]) + K2*d["cx"]**2
        ax[2].plot(d["cx"], n**1.5*d["cdE"] - cb, "o", ms=9, mfc="white", mec=colour[n], mew=2)
    for h in np.arange(0, a.xmax + 0.01, 0.5):
        for axh in ax:
            axh.axvline(h, color="#cccccc", lw=1, zorder=0)
    for axh in ax[1:]:
        axh.axhline(0, color="#888888", lw=1, zorder=0)
    ax[0].axhline(1/(4*math.sqrt(2*math.pi)), color="#bbbbbb", lw=1, ls=":", zorder=0)
    ax[0].set_ylim(-0.005, 0.13)
    ax[0].set_ylabel("$\\sqrt{n}$ (L $-$ L(1/2))")
    ax[0].set_title("(a)  L in the arch scaling; dashed: the limit A(x) = (1/4 $-$ 4 d(x)$^2$)/$\\sqrt{2\\pi}$, "
                    "d = distance to Z/2 + 1/4;  dotted: 1/(4$\\sqrt{2\\pi}$)", loc="left")
    ax[1].set_ylabel("n$^{1.5}$(L $-$ L(1/2)) $-$ n A(x) $-$ 2$\\sqrt{2/\\pi}$ x$^2$")
    ax[1].set_title("(b)  next order for L", loc="left")
    ax[2].set_ylabel("the same for E")
    ax[2].set_title("(c)  next order for E (circles: cusps)", loc="left")
    ax[2].set_xlabel("x = n(p $-$ 1/2)")
    ax[2].set_xlim(0, a.xmax)
    lo = min(np.percentile(n**1.5*D[n]["dL"] - n*arch(D[n]["x"]) - K2*D[n]["x"]**2, 0.5) for n in ns)
    hi = max(np.percentile(n**1.5*D[n]["dE"] - n*arch(D[n]["x"]) - K2*D[n]["x"]**2, 99.5) for n in ns)
    pad = 0.08*(hi - lo)
    ax[1].set_ylim(lo - pad, hi + pad); ax[2].set_ylim(lo - pad, hi + pad)
    hand = [Line2D([], [], color=colour[n], lw=3, label=f"n = {n}") for n in ns] + \
           [Line2D([], [], color="black", lw=2.2, ls="--", label="A(x)")]
    ax[0].legend(handles=hand, loc="upper right", ncol=2, fontsize=18, framealpha=0.92)

    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, f"L_arches_x{a.xmax:g}_n{'-'.join(map(str, ns))}.png")
    fig.savefig(path, dpi=100, bbox_inches="tight")
    print(path)

if __name__ == "__main__":
    main()
