"""plotting/lower_bound/_data.py -- E and L near p = 1/2 on an exact-kink grid, shared by the plots here.

L(p) = n + 1/2 - 2 D(p) is the certificate envelope of the (★) reduction (RESEARCH_LOG section 4 and
the 2026-09-24 claude.ai entry).  Nothing is derived here: E from binom_core.E_at / E_half, D from
star_check.D_direct, D(1/2) from star_check.D_half_exact, switch points from star_check.switch_points,
cusps from cusps/nNNNNN.csv.  The grid in x = n(p - 1/2) is uniform plus every tie point and switch
point in (0, xmax], plus optional dense windows of half-width `local`/n around each switch point, so
every kink of E and L is a grid point.  Descriptive, double precision.
"""
import csv, math, os, sys
from multiprocessing import Pool
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
import binom_core as core
from star_check import switch_points, D_direct, D_half_exact

OUT = os.path.join(ROOT, "lower_bound_proof_plots")
K2 = 2.0*math.sqrt(2.0/math.pi)                   # smooth term of the troughs: n^1.5 (L - L(1/2)) ~ K2 x^2

def arch(x):
    """Limit of sqrt(n) (L - L(1/2)) and of sqrt(n) (E - E(1/2)) at fixed x (RESEARCH_LOG 2026-10-03):
    (1/4 - 4 d^2)/sqrt(2 pi), d = distance from x to the nearest point of Z/2 + 1/4."""
    x = np.asarray(x, float)
    d = np.abs((x - 0.25) - 0.5*np.round((x - 0.25)/0.5))
    return (0.25 - 4.0*d*d)/math.sqrt(2.0*math.pi)

def trough_limit(m, n):
    """Limit of n^1.5 (L - L(1/2)) at the m-th switch point (RESEARCH_LOG 2026-09-24, result 4)."""
    m = np.asarray(m, float)
    lat = np.where(m % 2 == 1, np.where(n % 2 == 0, -0.25, 0.25), 0.0)
    return (m**2 + lat)/math.sqrt(2.0*math.pi)

def _tie_p(n, xmax):
    l = core.lnC_arr(n)
    i, j = np.triu_indices(n + 1, 1)
    keep = i + j > n
    i, j = i[keep], j[keep]
    p = 1.0/(1.0 + np.exp(-(l[i] - l[j])/(j - i)))
    return p[(p > 0.5) & (n*(p - 0.5) <= xmax)]

def _cusps(n, xmax):
    with open(os.path.join(ROOT, "cusps", f"n{n:05d}.csv")) as fh:
        rows = [(float(r["pstar"]), float(r["E"])) for r in csv.DictReader(fh)]
    p, E = np.array(rows).T
    k = n*(p - 0.5) <= xmax
    return p[k], E[k]

def _chunk(args):
    n, ps = args
    return [core.E_at(n, v) for v in ps], [D_direct(n, v) for v in ps]

def compute(ns, xmax, grid=3000, local=None, workers=8):
    """{n: dict} with
        x            grid in x = n(p - 1/2)
        dE, dL       E - E(1/2) and L - L(1/2) = 2 (D(1/2) - D) on the grid (unscaled)
        cx, cdE      cusps of E in range: x and E - E(1/2)
        m, xm, dLm   switch points in range: index, x, L - L(1/2)"""
    plan = {}
    for n in ns:
        m, _, _, ps, _ = switch_points(n)
        k = n*(ps - 0.5) <= xmax
        xs = [np.linspace(0, xmax, grid + 1)[1:]]
        if local:
            xs += [n*(v - 0.5) + np.linspace(-local, local, 401)/n for v in ps[k]]
        x = np.concatenate(xs)
        x = x[(x > 0) & (x <= xmax)]
        p = np.unique(np.concatenate([0.5 + x/n, _tie_p(n, xmax), ps[k]]))
        plan[n] = (p, m[k], ps[k])
    jobs = [(n, c) for n in ns for c in np.array_split(plan[n][0], max(1, len(plan[n][0])//200))]
    with Pool(workers) as pool:
        res = pool.map(_chunk, jobs)
    raw = {n: ([], []) for n in ns}
    for (n, _), (e, d) in zip(jobs, res):
        raw[n][0].extend(e); raw[n][1].extend(d)
    out = {}
    for n in ns:
        p, m, ps = plan[n]
        x = n*(p - 0.5)
        dE = np.array(raw[n][0]) - core.E_half(n)
        dL = 2.0*(D_half_exact(n) - np.array(raw[n][1]))
        cp, cE = _cusps(n, xmax)
        xm = n*(ps - 0.5)
        out[n] = dict(x=x, dE=dE, dL=dL, cx=n*(cp - 0.5), cdE=cE - core.E_half(n),
                      m=m, xm=xm, dLm=np.interp(xm, x, dL))
        print(f"n={n}: {len(p)} points")
    return out
