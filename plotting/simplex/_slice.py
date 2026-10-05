"""plotting/simplex/_slice.py -- exact 3D slices of the n-simplex through the touch point p = 1/2.

For n >= 4 the simplex cannot be drawn whole, and a projection would replace B by its shadow, which
covers the curve.  So we cut instead.  With g(t) = f at p = 1/2 + t, the symmetry x_k <-> x_{n-k}
makes every even derivative of g at 0 a symmetric vector and every odd one antisymmetric.  The flat

    x = x* + A u,   x* = f_{1/2},   A = [tangent g'(0) | symmetric | -g''(0)]  (orthonormal, sum-zero)

therefore contains g(0), g'(0) and g''(0): the curve lies in it to second order and leaves it only
through g'''.  For n = 4 and 5 the symmetric sum-zero vectors form a 2D space, so the flat is
span(symmetric) + tangent exactly; u = (along the curve, the other symmetric direction, away from the
centre).  Every region is still exact in the slice: a half-space <a, x> <= t pulls back to
<A^T a, u> <= t - <a, x*>.  The curve is drawn as its orthogonal projection A^T (f - x*), and
off(p) = |f - x* - A A^T (f - x*)| / |f - x*| says how far it has left the slice.
"""
import math
import numpy as np
from scipy.spatial import HalfspaceIntersection

from _simplex import binomial
from _simplex3d import hull

class Slice:
    def __init__(self, n):
        self.n = n
        k = np.arange(n + 1)
        self.xs = binomial(n, [0.5])[0]
        t = np.array([math.comb(n, j) for j in k]) * (2*k - n)            # g'(0), up to a factor
        g2 = np.array([math.comb(n, j) for j in k]) * ((2*k - n)**2 - n)   # g''(0), up to a factor
        t = t/np.linalg.norm(t)
        z = -g2/np.linalg.norm(g2)                                        # away from the centre
        sym = [np.eye(n + 1)[j] + np.eye(n + 1)[n - j] for j in range(n//2 + 1)]
        sym = [v - v.sum()/(n + 1) for v in sym]                          # sum-zero symmetric vectors
        y = None
        for v in sym:                                                     # first one independent of z
            w = v - (v @ z)*z
            if np.linalg.norm(w) > 1e-9: y = w/np.linalg.norm(w); break
        self.A = np.c_[t, y, z]
        self.centre_u = self.u(np.full(n + 1, 1.0/(n + 1)))

    def u(self, x):
        return (np.asarray(x) - self.xs) @ self.A

    def off(self, x):
        d = np.asarray(x) - self.xs
        r = d - (d @ self.A) @ self.A.T
        return np.linalg.norm(r, axis=-1)/np.maximum(np.linalg.norm(d, axis=-1), 1e-300)

    def polytope(self, halfspaces, interior_u=None):
        """(V, faces) in u-coordinates of {<a, x> <= t} intersected with the simplex and the flat."""
        A, xs = self.A, self.xs
        H = [np.r_[A.T @ a, a @ xs - t] for a, t in halfspaces]
        H += [np.r_[-A[i], -xs[i]] for i in range(self.n + 1)]
        ip = self.centre_u if interior_u is None else interior_u
        hsi = HalfspaceIntersection(np.array(H), ip)
        return hull(np.unique(np.round(hsi.intersections, 12), axis=0))

    def plane_in_slice(self, a, t):
        """{<a, x> = t} within the simplex slice, as a polygon in u (or None if it misses)."""
        V, F = self.polytope([])
        b, c = self.A.T @ a, t - a @ self.xs
        if np.linalg.norm(b) < 1e-12: return None
        s = V @ b - c
        if s.max() <= 1e-12 or s.min() >= -1e-12: return None
        pts = []
        for f in F:                                                       # cut every edge
            for i, j in zip(f, f[1:] + f[:1]):
                if s[i]*s[j] < 0: pts.append(V[i] + (V[j] - V[i])*s[i]/(s[i] - s[j]))
                elif abs(s[i]) < 1e-12: pts.append(V[i])
        P = np.unique(np.round(pts, 12), axis=0)
        cen = P.mean(0); nrm = b/np.linalg.norm(b)
        e1 = P[0] - cen; e1 /= np.linalg.norm(e1); e2 = np.cross(nrm, e1)
        return P[np.argsort(np.arctan2((P - cen) @ e2, (P - cen) @ e1))]


def section(V, faces, axis=1, value=0.0):
    """Cut a 3D polytope (V, faces) by the plane u[axis] = value: the section polygon, as its two
    remaining coordinates in order (or None)."""
    s = V[:, axis] - value
    pts = []
    for f in faces:
        for i, j in zip(f, f[1:] + f[:1]):
            if s[i]*s[j] < 0: pts.append(V[i] + (V[j] - V[i])*s[i]/(s[i] - s[j]))
            elif abs(s[i]) < 1e-12: pts.append(V[i])
    if len(pts) < 3: return None
    keep = [k for k in range(3) if k != axis]
    P = np.unique(np.round(np.array(pts)[:, keep], 12), axis=0)
    c = P.mean(0)
    return P[np.argsort(np.arctan2(P[:, 1] - c[1], P[:, 0] - c[0]))]
