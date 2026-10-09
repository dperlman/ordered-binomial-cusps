#!/usr/bin/env python3
"""The cusp curve in pair space, and F3 < 0 as the cusps off it: checks behind the 2026-10-09 entry.

Pair space: each tie point (i, j) at grid position x = (i+j+1)/(2(n+1)) (~ p*, RESEARCH_LOG fact 12)
and scaled width w = (j-i)/sqrt(n); u = S_-/kappa is where E's slope changes sign relative to the
kink (cusp <=> -1 < u < 0; u <= -1: E falls through the tie point; u >= 0: it rises).  Prints:
  1. the main cusp curve: median w of the cusps with w < 2 at several x, for several n, from
     cusps/cusps_all.csv, its last x, and the share of cusps on it;
  2. the cusps off the curve (w >= 2) against the cusps with F3 < 0, over every n <= 5000;
  3. along each column of fixed i+j (narrow tie points, w < 2), from full tie tables
     (obd_core.tie_table): E falls through the narrowest tie points and switches once to rising,
     and the column's cusp is the first tie point that is not falling.

    .venv/bin/python cusp_curve.py            # about a minute (tie tables for n = 500..4000)
"""
from multiprocessing import Pool

import numpy as np
import pyarrow.csv as pc

import obd_core

CSV = "cusps/cusps_all.csv"


def catalogue():
    t = pc.read_csv(CSV, convert_options=pc.ConvertOptions(include_columns=["n", "i", "j", "F3", "F3_sign"]))
    n, i, j, f3 = (t[c].to_numpy().astype(float) for c in ("n", "i", "j", "F3"))
    neg = t["F3_sign"].to_numpy(zero_copy_only=False) == "-"
    x = (i + j + 1) / (2 * (n + 1))
    w = (j - i) / np.sqrt(n)

    print("1. main cusp curve: median width/sqrt(n) of the cusps narrower than 2 sqrt(n), at grid position x")
    xs = (0.505, 0.55, 0.60, 0.63, 0.645)
    print("     n   " + "  ".join(f"{b:.3f}" for b in xs) + "   last x    share of the cusps of n")
    for N in (1000, 2000, 3000, 4000, 5000):
        s = (n == N) & (w < 2)
        cells = []
        for b in xs:
            k = s & (np.abs(x - b) < 0.0025)
            cells.append(f"{np.median(w[k]):.3f}" if k.any() else "  -  ")
        print(f"  {N:5d}  " + "  ".join(cells) + f"   {x[s].max():.4f}   {s.sum() / (n == N).sum():.1%}")

    print("\n2. the cusps off the curve against F3 < 0, every n = 3..5000")
    wide = w >= 2
    print(f"   F3 < 0: {neg.sum():,} cusps; width >= 2 sqrt(n): {wide.sum():,}; both: {(neg & wide).sum():,}")
    for q in np.flatnonzero(wide != neg):
        print(f"   exception: n = {int(n[q])}, pair ({int(i[q])}, {int(j[q])}), width/sqrt(n) = {w[q]:.2f}, F3 = {f3[q]:+.4g}")
    for lo, hi in [(3, 99), (100, 999), (1000, 5000)]:
        s = (n >= lo) & (n <= hi)
        print(f"   n = {lo}..{hi}: F3 >= 0 widths up to {w[s & ~neg].max():.3f} sqrt(n); "
              f"F3 < 0 widths from {w[s & neg].min():.3f} sqrt(n)")


def columns(pool):
    print("\n3. columns of fixed i+j, narrow tie points (width < 2 sqrt(n)), grid position < 0.66, from tie_table")
    for n in (500, 1000, 2000, 4000):
        t = obd_core.tie_table(n, workers=8, pool=pool)
        s = slice(1, None)                                       # row 0 is the axis p = 1/2
        i, j = t["i"][s], t["j"][s]
        m = j - i
        sm, cusp = t["S_minus"][s], t["is_cusp"][s]
        lu = np.log(np.abs(sm) + 1e-300) - np.log(m) - t["ln_fi"][s]
        u = np.sign(sm) * np.exp(np.clip(lu, -700, 700))     # u = S_-/kappa in logs: no overflow
        col = i + j
        sel = (m < 2 * np.sqrt(n)) & ((col + 1) / (2 * (n + 1)) < 0.66)
        o = np.lexsort((m[sel], col[sel]))
        c, uu, cu = col[sel][o], u[sel][o], cusp[sel][o]
        fall = uu <= -1
        starts = np.r_[0, np.flatnonzero(np.diff(c)) + 1]
        ends = np.r_[starts[1:], c.size]
        shape = at_switch = 0
        zero_x, zero_all_rising = [], 0
        per_col = np.zeros(4, int)                               # columns with 0, 1, 2, 3+ narrow cusps
        for a, b in zip(starts, ends):
            f = fall[a:b]
            sw = np.flatnonzero(f[:-1] != f[1:])
            shape += (sw.size == 0) or (sw.size == 1 and f[0])
            cc = np.flatnonzero(cu[a:b])
            per_col[min(cc.size, 3)] += 1
            if cc.size == 0:
                zero_x.append((c[a] + 1) / (2 * (n + 1)))
                zero_all_rising += not f.any()
            nf = np.flatnonzero(~f)
            # the cusps are exactly the first tie points after the switch, consecutively
            at_switch += cc.size > 0 and nf.size > 0 and np.array_equal(cc, np.arange(nf[0], nf[0] + cc.size))
        print(f"   n = {n}: {starts.size} columns; falling then rising, at most one switch: {shape}; "
              f"columns with 0/1/2/3+ narrow cusps: {'/'.join(map(str, per_col))}; "
              f"cusps exactly the first tie points after the switch: {at_switch} of {per_col[1:].sum()}; "
              f"narrow cusps {int(cu.sum())} of {int(cusp.sum())}")
        if zero_x:
            print(f"      columns with no narrow cusp: grid positions {min(zero_x):.4f}..{max(zero_x):.4f}, "
                  f"{zero_all_rising} of {len(zero_x)} with E rising through every narrow tie point (no switch)")


if __name__ == "__main__":
    catalogue()
    with Pool(8) as pool:
        columns(pool)
