"""
E at the first, second and third cusp of each n, against n.

    .venv/bin/python plotting/first_cusps.py [--csv cusps/cusps_all.csv] [--k 3]

"First" means lowest p* above 1/2.  Three panels, because raw E hides everything:

  (a) E itself.  E ~ n/2, so the three curves are visually identical -- at n=1000 all three read
      975.2624.  The panel is here to show that, not to be read.
  (b) E - E(1/2), log-log.  Now they separate into parallel lines of slope -3/2.
  (c) (E - E(1/2)) * n^(3/2).  Flat, which is the point: each series has a constant.

Only n with at least three cusps are drawn (all n except 4, 5, 7, 8), and only n <= 3000, where
every n is covered.  The n = 4000..8000 Parquet dumps are deliberately excluded: five scattered
values of n would read as a trend they cannot support.
"""
import argparse, csv, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import binom_core as core

COL = ["#2b5d8a", "#e08b2d", "#1b7a4b", "#8a4fbf", "#c1272d"]

def load(path, k):
    byn = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            if r['pstar'] == '': continue
            byn.setdefault(int(r['n']), []).append((float(r['pstar']), float(r['E'])))
    ns = np.array([n for n in sorted(byn) if len(byn[n]) >= k])
    E = np.array([[e for _, e in sorted(byn[n])[:k]] for n in ns])       # ordered by p*
    Eh = np.array([core.E_half(n) for n in ns])
    return ns, E, Eh

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="cusps/cusps_all.csv")
    ap.add_argument("--out", default="plots"); ap.add_argument("--k", type=int, default=3)
    ap.add_argument("--width", type=int, default=6000); ap.add_argument("--height", type=int, default=3000)
    a = ap.parse_args()
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ns, E, Eh = load(a.csv, a.k)
    d = E - Eh[:, None]
    scaled = d * (ns[:, None]**1.5)
    ordin = ["1st", "2nd", "3rd", "4th", "5th"]
    dpi = 100
    fig, axes = plt.subplots(3, 1, figsize=(a.width/dpi, a.height/dpi), dpi=dpi,
                             gridspec_kw=dict(hspace=0.30))

    ax = axes[0]
    ax.plot(ns, Eh, lw=3.5, color="0.45", ls="--", label="$E(n,1/2)$")
    for t in range(a.k):
        ax.plot(ns, E[:, t], lw=2.2, color=COL[t], label=f"E at the {ordin[t]} cusp")
    ax.set_xscale("log"); ax.set_ylabel("$E(n,p^*)$", fontsize=32, labelpad=18)
    ax.set_title("(a)  E itself — the three curves and $E(n,1/2)$ are indistinguishable "
                 "(all four read 975.2624 at n=1000)", fontsize=34, pad=18, loc="left")
    ax.legend(fontsize=26, loc="upper left", framealpha=0.92)

    ax = axes[1]
    for t in range(a.k):
        for par, ls, lab in ((0, "-", "even n"), (1, ":", "odd n")):
            m = ns % 2 == par
            ax.plot(ns[m], d[m, t], lw=2.0, ls=ls, color=COL[t],
                    label=f"{ordin[t]} cusp, {lab}")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_ylabel(r"$E(n,p^*) - E(n,1/2)$", fontsize=32, labelpad=18)
    ax.set_title(r"(b)  the difference — parallel lines of slope $-3/2$; each splits by the parity of $n$",
                 fontsize=34, pad=18, loc="left")
    ax.legend(fontsize=22, loc="lower left", ncol=a.k, framealpha=0.92)

    ax = axes[2]
    txt = []
    for t in range(a.k):
        for par, ls, lab in ((0, "-", "even"), (1, ":", "odd")):
            m = (ns % 2 == par) & (ns >= 500)
            c = scaled[m, t]
            ax.plot(ns[ns % 2 == par], scaled[ns % 2 == par, t], lw=2.0, ls=ls, color=COL[t])
            txt.append(f"{ordin[t]} {lab}: {c.mean():.4f}")
    ax.set_xscale("log")
    ax.set_ylabel(r"$(E-E(1/2))\,\cdot\,n^{3/2}$", fontsize=32, labelpad=18)
    ax.set_xlabel("$n$", fontsize=36, labelpad=18)
    ax.set_title("(c)  scaled by $n^{3/2}$ — flat, so each series is a constant.   "
                 + ";   ".join(txt), fontsize=30, pad=18, loc="left")
    for x in axes:
        x.tick_params(labelsize=24, length=10, width=2); x.grid(True, alpha=0.25, lw=1.1)
        x.set_xlim(ns.min()*0.9, ns.max()*1.05)
    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, f"first_cusps_E.png")
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    print(f"{path}  ({os.path.getsize(path)/1e6:.1f} MB);  n = {ns.min()}..{ns.max()}, "
          f"{len(ns)} values ({', '.join(str(x) for x in [4,5,7,8])} have fewer than {a.k} cusps)")
    for t in range(a.k):
        for par, lab in ((0, "even"), (1, "odd ")):
            m = (ns % 2 == par) & (ns >= 500)
            print(f"  {ordin[t]} cusp, {lab} n:  (E-E(1/2))*n^1.5 = {scaled[m,t].mean():.5f} "
                  f"+- {scaled[m,t].std():.1e}")

if __name__ == "__main__":
    main()
