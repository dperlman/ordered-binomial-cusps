"""plotting/simplex/_simplex.py -- shared geometry for the simplex figures.

The masses x = (f_p(0), ..., f_p(n)) live in the probability simplex.  E(x) = sum_k k * sorted(x)_k
is the largest of the (n+1)! linear functionals <sigma, x> (rearrangement inequality): the support
function of the permutohedron.  It is linear on each CHAMBER (a fixed ordering of the coordinates);
the chamber walls x_i = x_j are crossed at the tie points.  A certificate W_c (RESEARCH_LOG section 4)
is the linear functional <a_c, x> with a_c[k] = n + 1/2 - 2|k - c|, c in Z/2 + 1/4, and L = max_c W_c.

Polygons are lists of barycentric points (length n+1 arrays) clipped by half-spaces <a, x> <= t,
so every region drawn is exact.  Output directory: simplex_figures/ (gitignored).
"""
import itertools, math, os
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "simplex_figures")

# ---- the triangle (n = 2): e0 bottom left, e1 top, e2 bottom right
TRI = np.array([[0.0, 0.0], [0.5, math.sqrt(3)/2], [1.0, 0.0]])

def xy(x):
    """Barycentric (..., 3) -> plane (..., 2) for the triangle."""
    return np.asarray(x) @ TRI

def clip(poly, a, t):
    """Sutherland-Hodgman: the part of the convex polygon `poly` with <a, x> <= t."""
    out = []
    for k in range(len(poly)):
        P, Q = poly[k], poly[(k + 1) % len(poly)]
        fp, fq = a @ P - t, a @ Q - t
        if fp <= 0: out.append(P)
        if fp*fq < 0: out.append(P + (Q - P)*fp/(fp - fq))
    return out

def simplex(n):
    return [np.eye(n + 1)[i] for i in range(n + 1)]

def region(n, halfspaces):
    """Intersection of the simplex with every <a, x> <= t in `halfspaces`."""
    poly = simplex(n)
    for a, t in halfspaces:
        poly = clip(poly, np.asarray(a, float), t)
        if not poly: break
    return poly

def chamber(order):
    """Chamber x[order[0]] >= x[order[1]] >= ... as a polygon."""
    n = len(order) - 1
    hs = []
    for hi, lo in zip(order[:-1], order[1:]):
        a = np.zeros(n + 1); a[lo] = 1; a[hi] = -1          # x_lo - x_hi <= 0
        hs.append((a, 0.0))
    return region(n, hs)

def perms(n):
    return [np.array(s, float) for s in itertools.permutations(range(n + 1))]

def E_sublevel(n, t):
    """{x : E(x) <= t}: one facet <sigma, x> <= t per permutation sigma."""
    return region(n, [(s, t) for s in perms(n)])

def certificates(n):
    """{c: a_c} for c = -1/4 .. n + 1/4 (the same family as L_small_n.py)."""
    cs = np.arange(-1, 2*n + 1)/2.0 + 0.25
    k = np.arange(n + 1)
    return {float(c): n + 0.5 - 2.0*np.abs(k - c) for c in cs}

def L_sublevel(n, t):
    return region(n, [(a, t) for a in certificates(n).values()])

def binomial(n, p):
    """Mass vectors f_p, shape (len(p), n+1)."""
    p = np.asarray(p, float)[:, None]
    k = np.arange(n + 1)
    return np.array([math.comb(n, j) for j in k]) * p**k * (1 - p)**(n - k)

def E_of(x):
    return np.sort(x, axis=-1) @ np.arange(x.shape[-1])

def line_in_simplex(n, a, t):
    """The segment of {<a, x> = t} inside the simplex (two barycentric endpoints, or [])."""
    poly = simplex(n)
    pts = []
    for k in range(len(poly)):
        P, Q = poly[k], poly[(k + 1) % len(poly)]
        fp, fq = a @ P - t, a @ Q - t
        if fp == 0: pts.append(P)
        elif fp*fq < 0: pts.append(P + (Q - P)*fp/(fp - fq))
    return pts[:2]
