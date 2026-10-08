#!/usr/bin/env python3
"""Cusps sit next to the fractions k/(2(n+1)): checks behind the 2026-10-08 RESEARCH_LOG entry.

Reads only cusps/cusps_all.csv (n = 3..5000, p* > 1/2).  Prints:
  1. the identity  logit p* = mean of logit(t/(n+1)), t = i+1..j  on a sample of cusps;
  2. the offset of p* above the grid point (i+j+1)/(2(n+1)), against its second-order (Jensen) form;
  3. FCW_r(p) = the first n with a cusp within r of p: its median times sqrt(r);
  4. the spikes of FCW beside fractions a/b: height against the Farey gap, and left/right asymmetry;
  5. which cusp first comes within r just left / right of a few a/b.
FCW uses every cusp of n = 3..5000 plus, for every n >= 2, the axis cusp p = 1/2 (always a cusp;
n = 2 has no other) and the mirror images 1 - p*, as the OBD cusp-proximity plots do.

    .venv/bin/python farey_cusps.py          # about a minute
"""
import math

import numpy as np
import pyarrow.csv as pc

CSV = "cusps/cusps_all.csv"


def first_n_within_r(cn, cp, grid, r):
    """Smallest n with a cusp within r of each grid p (0: none).  cn, cp sorted by n."""
    lo = np.searchsorted(grid, cp - r, "left")
    hi = np.searchsorted(grid, cp + r, "right")
    keep = hi > lo
    lo, hi, nn = lo[keep], hi[keep], cn[keep]
    best = np.full(grid.size, np.iinfo(np.int64).max)
    for k in range(int((hi - lo).max())):
        m = lo + k < hi
        np.minimum.at(best, lo[m] + k, nn[m])
    best[best == np.iinfo(np.int64).max] = 0
    return best


def main():
    t = pc.read_csv(CSV, convert_options=pc.ConvertOptions(include_columns=["n", "i", "j", "pstar"]))
    n, i, j, p = (t[c].to_numpy() for c in ("n", "i", "j", "pstar"))
    nf, i_f, j_f = n.astype(float), i.astype(float), j.astype(float)
    m = j_f - i_f
    print(f"{n.size} cusps, n = {n.min()}..{n.max()}\n")

    # 1. the identity
    rng = np.random.default_rng(0)
    worst = 0.0
    for k in rng.integers(0, n.size, 20000):
        tt = np.arange(i[k] + 1, j[k] + 1) / (n[k] + 1.0)
        lp = float(np.mean(np.log(tt) - np.log1p(-tt)))
        worst = max(worst, abs(1 / (1 + math.exp(-lp)) - p[k]))
    print(f"1. logit p* = mean logit(t/(n+1)): worst |difference| over 20,000 random cusps {worst:.1e}")

    # 2. offset above the grid point, in grid steps 1/(2(n+1))
    mid = (i_f + j_f + 1) / (2 * (nf + 1))
    off = 2 * (nf + 1) * (p - mid)
    pred = (2 * mid - 1) / (mid * (1 - mid)) * (m ** 2 - 1) / (12 * (nf + 1))
    s = n >= 100
    rel = np.abs(off - pred) / np.maximum(pred, 1e-300)
    print(f"2. offset >= 0 for {np.mean(off[s] >= -1e-9):.4%} of cusps with n >= 100; "
          f"second-order form: median relative error {np.median(rel[s & (m > 2)]):.2%}, "
          f"99th percentile {np.quantile(rel[s & (m > 2)], .99):.2%}")
    q = np.quantile(off[s], [.5, .9, .99, 1.0])
    print("   offset in grid steps: median %.3f, 90%% %.3f, 99%% %.3f, max %.3f" % tuple(q))
    print("   width/sqrt(n): median %.2f, 99%% %.2f, max %.2f" % tuple(np.quantile(m[s] / np.sqrt(nf[s]), [.5, .99, 1])))

    # FCW data: cusps, mirrors, and the axis for every n >= 2
    axis_n = np.arange(2, n.max() + 1)
    cn = np.concatenate([n, n, axis_n]).astype(np.int64)
    cp = np.concatenate([p, 1 - p, np.full(axis_n.size, 0.5)])
    o = np.lexsort((cp, cn))
    cn, cp = cn[o], cp[o]

    # 3. median FCW * sqrt(r) in the band
    print("\n3. median of FCW_r(p) * sqrt(r) over 0.51 < p < 0.65:")
    for r in (1e-3, 1e-4, 1e-5, 1e-6):
        g = np.linspace(0.51, 0.65, 200001)
        f = first_n_within_r(cn, cp, g, r)
        print(f"   r = {r:g}: {np.median(f[f > 0]) * math.sqrt(r):.3f}")

    # 4. spikes beside a/b, r = 1e-5, grid spacing r/10
    r = 1e-5
    g = np.linspace(0.5, 0.657, int(round(0.157 / (r / 10))) + 1)
    f = first_n_within_r(cn, cp, g, r).astype(float)
    f[f == 0] = n.max() + 1
    rows = []
    for b in range(5, 41):
        for a in range(math.ceil(0.505 * b), math.floor(0.648 * b) + 1):
            if math.gcd(a, b) != 1:
                continue
            x = a / b
            local = np.median(f[(np.abs(g - x) > 0.003) & (np.abs(g - x) < 0.006)])
            left = f[(g < x - r) & (g > x - 3e-4)].max() / local
            right = f[(g > x + r) & (g < x + 3e-4)].max() / local
            gap = (2 if b % 2 == 0 else 1) / b       # Farey gap around a/b on the grid, in grid steps
            rows.append((b, gap, left, right))
    rows = np.array(rows)
    b_, gap, left, right = rows.T
    print("\n4. spikes at r = 1e-5 (height = max FCW within 3e-4 beside a/b, over the local median):")
    print(f"   corr(log left height, log Farey gap) = {np.corrcoef(np.log(left), np.log(gap))[0, 1]:.2f}; "
          f"against 1/b alone {np.corrcoef(np.log(left), -np.log(b_))[0, 1]:.2f}")
    for lo, hi in [(5, 10), (11, 20), (21, 28), (29, 40)]:
        s = (b_ >= lo) & (b_ <= hi)
        print(f"   b = {lo}-{hi}: median left height {np.median(left[s]):.1f}x, right/left {np.median(right[s] / left[s]):.2f}")

    # 5. which cusp first comes within r just beside a/b
    print("\n5. the cusp that first comes within r = 1e-5 of p, for 40 p from r to r + 2e-5 beside a/b:")
    for a, b in [(3, 5), (4, 7), (5, 8), (7, 12), (9, 16)]:
        x = a / b
        for side, sg in (("left ", -1), ("right", 1)):
            first, on_grid = [], 0
            for pp in x + sg * np.linspace(r + 1e-6, r + 2e-5, 40):
                near = np.flatnonzero(np.abs(p - pp) <= r)
                k = near[np.argmin(n[near])]
                first.append(n[k])
                # is its grid point (i+j+1)/(2(n+1)) the fraction a/b itself?
                on_grid += (int(i[k] + j[k] + 1) * b == 2 * (int(n[k]) + 1) * a)
            print(f"   {a}/{b} {side}: median first n {int(np.median(first)):5d}; "
                  f"grid point = {a}/{b} in {on_grid}/40")


if __name__ == "__main__":
    main()
