"""plotting/_style.py -- shared figure defaults for every plot in this project.

DEFAULT RESOLUTION IS 6000x3000 (changed 2026-09-22 from 10000x6000, which was far larger than
anything needed and slow to open).  Pass --width/--height to override.

ALIASING.  These plots routinely put ~5000 points on the x axis, which at 6000 px is about 1.1 px
per point.  Drawing a 1 px marker at 1.1 px spacing produces moire banding: the marker lattice and
the pixel lattice beat against each other.  The fix is to make markers SLIGHTLY WIDER than the
point spacing so neighbours overlap and blend, keep alpha below 1 so the overlap averages rather
than replaces, leave antialiasing on, and NEVER connect dense points with a line (a polyline
through 5000 points at 1 px spacing is pure aliasing).  marker_size() below does the first part.
"""
import numpy as np

DPI = 100
DEFAULT_W, DEFAULT_H = 6000, 3000

def add_size_args(ap, w=DEFAULT_W, h=DEFAULT_H):
    ap.add_argument("--width", type=int, default=w)
    ap.add_argument("--height", type=int, default=h)
    return ap

def figsize(a):
    return (a.width/DPI, a.height/DPI)

def marker_size(n_points, axes_frac=0.88, width_px=DEFAULT_W, overlap=3.2):
    """Marker area in points^2 such that markers are ~`overlap` x the point spacing.

    At 6000 px and 5000 points the spacing is ~1.1 px, so a ~3.5 px marker overlaps its neighbours
    and the series reads as a continuous band instead of a moire pattern.  Clamped so that sparse
    series still get visible dots and dense ones do not turn into a solid slab.
    """
    spacing_px = max(width_px*axes_frac/max(n_points, 1), 0.05)
    d_px = float(np.clip(spacing_px*overlap, 3.0, 26.0))
    return (d_px/(DPI/72.0))**2

def save(fig, path, dpi=DPI):
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    return path
