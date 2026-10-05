"""plotting/simplex/_simplex3d.py -- the n = 3 simplex (a tetrahedron) in R^3, and exact polytopes in it.

Embedding: e0 = (-1, 0, -s), e3 = (1, 0, -s), e1 = (0, -1, s), e2 = (0, 1, s), s = 1/sqrt(2): a regular
tetrahedron with edge 2.  The symmetry p <-> 1-p (x_k <-> x_{3-k}) is the half-turn about the z axis,
and the binomial curve runs from e0 (p = 0) to e3 (p = 1), meeting the z axis at p = 1/2.
A polytope is given by half-spaces <a, x> <= t in barycentric coordinates (plus x >= 0); its vertices
come from scipy's HalfspaceIntersection and its faces from ConvexHull, coplanar triangles merged.
"""
import math
import numpy as np
from scipy.spatial import ConvexHull, HalfspaceIntersection

S = 1/math.sqrt(2)
TET = np.array([[-1.0, 0.0, -S], [0.0, -1.0, S], [0.0, 1.0, S], [1.0, 0.0, -S]])
_M = np.linalg.inv(np.vstack([TET.T, np.ones(4)]))          # x = _M @ [X, 1]

def xyz(x):
    """Barycentric (..., 4) -> R^3 (..., 3)."""
    return np.asarray(x) @ TET

def polytope(halfspaces, interior):
    """Vertices (k, 3) and faces (list of index loops) of {<a, x> <= t for all} within the simplex.
    `interior` is a barycentric point strictly inside."""
    rows = [(np.asarray(a, float), float(t)) for a, t in halfspaces]
    rows += [(-np.eye(4)[k], 0.0) for k in range(4)]
    H = np.array([np.r_[a @ _M[:, :3], a @ _M[:, 3] - t] for a, t in rows])
    hsi = HalfspaceIntersection(H, xyz(interior))
    return hull(np.unique(np.round(hsi.intersections, 12), axis=0))   # a vertex on >3 planes repeats

def hull(V):
    """Vertices (k, 3) -> (V, faces): convex hull with coplanar triangles merged into polygons."""
    hull = ConvexHull(V)
    groups = {}
    for simp, eq in zip(hull.simplices, hull.equations):
        key = tuple(np.round(eq, 7))
        groups.setdefault(key, set()).update(simp.tolist())
    faces = []
    for key, idx in groups.items():
        idx = list(idx)
        P = V[idx]; c = P.mean(0); nrm = np.array(key[:3])
        u = P[0] - c; u /= np.linalg.norm(u); w = np.cross(nrm, u)
        ang = np.arctan2((P - c) @ w, (P - c) @ u)
        faces.append([idx[k] for k in np.argsort(ang)])
    return V, faces

def face_on_plane(V, faces, a, t):
    """The face of a polytope lying on <a, x> = t (as an R^3 polygon), or None."""
    for f in faces:
        x = np.hstack([V[f], np.ones((len(f), 1))]) @ _M.T
        if np.all(np.abs(x @ a - t) < 1e-9): return V[f]
    return None

def edges(V, faces):
    """Unique edges of the merged faces, as (2, 3) segments."""
    seen, out = set(), []
    for f in faces:
        for i, j in zip(f, f[1:] + f[:1]):
            k = (min(i, j), max(i, j))
            if k not in seen:
                seen.add(k); out.append(V[[i, j]])
    return out

def triangles(faces):
    """Fan-triangulate merged faces: (i, j, k) index triples."""
    return [(f[0], f[m], f[m + 1]) for f in faces for m in range(1, len(f) - 1)]
