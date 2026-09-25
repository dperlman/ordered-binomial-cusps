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


# ---------------------------------------------------------------------------------------------
# ONE PIXEL COLUMN PER n  (the default for every plot with n on a LINEAR x axis, from 2026-09-23)
#
# Every n gets exactly k whole pixel columns, so every per-n mark is exactly k px wide and none is
# split across two columns.  With an integer spacing there is nothing for moire to beat against,
# so marker_size() and its overlap/alpha workaround do not apply here.  How it is done:
#   * the data area is placed in exact pixels (fixed pixel margins; the figure width FOLLOWS the
#     n-range: margins + N*k).  NEVER bbox_inches="tight" -- it re-crops after layout.
#   * xlim = (nmin-0.5, nmax+0.5), so n sits at the centre of its own column(s).
#   * the data are not drawn by matplotlib at all: they are painted into an RGBA array, one
#     column block per n, and placed with fig.figimage at an integer pixel offset -- no
#     resampling anywhere.  Rows come from the axes' own transData, so a log y axis just works.
#   * layering: a BACK axes (grid, ticks, labels; spines pushed 1 px outside the data area),
#     then the painted data, then a transparent FRONT axes for fitted curves, annotations and the
#     legend.  Draw overlays on g.ax / g.axes[i]; paint the data with g.points / g.density.
#   * k defaults to the largest integer with N*k <= DEFAULT_W, i.e. 1 px per n once N > 3000.
# Only for a linear n axis: --logx cannot have uniform columns, and uses marker_size() instead.
# View at 100%: any viewer that fits the image to the window resamples it and brings moire back.
# ---------------------------------------------------------------------------------------------
N_MARGINS = dict(left=300, right=90, bottom=190, top=170)   # px; room for the 28-44 pt text
SPINE_PX = 2

def _hex_rgb(c):
    import matplotlib.colors as mc
    return np.array(mc.to_rgb(c), np.float32)

class NGrid:
    def __init__(self, nmin, nmax, height=DEFAULT_H, k=None, nrows=1, gap=150,
                 margins=N_MARGINS, target_width=DEFAULT_W):
        import matplotlib.pyplot as plt
        self.nmin, self.nmax = int(nmin), int(nmax)
        self.N = self.nmax - self.nmin + 1
        self.k = int(k) if k else max(1, target_width // self.N)
        L, R, B, T = (margins[s] for s in ("left", "right", "bottom", "top"))
        self.pw = self.N*self.k
        self.ph = (height - B - T - gap*(nrows-1)) // nrows
        self.W, self.H = L + self.pw + R, B + T + nrows*self.ph + gap*(nrows-1)
        self.x0 = L
        self.y0 = [B + (nrows-1-r)*(self.ph + gap) for r in range(nrows)]   # row 0 on top
        self.fig = plt.figure(figsize=(self.W/DPI, self.H/DPI), dpi=DPI)
        off = (SPINE_PX/2 + 1)*72/DPI
        self.back, self.axes, self.canvas, self._ylim = [], [], [], [None]*nrows
        for r in range(nrows):
            rect = [L/self.W, self.y0[r]/self.H, self.pw/self.W, self.ph/self.H]
            b = self.fig.add_axes(rect, zorder=0)
            b.set_xlim(self.nmin - 0.5, self.nmax + 0.5)
            for s in b.spines.values():
                s.set_position(("outward", off)); s.set_linewidth(SPINE_PX*72/DPI)
            b.tick_params(direction="out")
            f = self.fig.add_axes(rect, sharex=b, sharey=b, zorder=2)
            f.patch.set_visible(False); f.axis("off")
            self.back.append(b); self.axes.append(f)
            self.canvas.append(np.zeros((self.ph, self.pw, 4), np.float32))   # premultiplied
        self.ax = self.axes[0]

    def _cols_rows(self, n, y, panel):
        n = np.asarray(n); y = np.asarray(y, float)
        if not np.issubdtype(n.dtype, np.integer) or n.min() < self.nmin or n.max() > self.nmax:
            raise ValueError("n must be integers inside the grid's range")
        b = self.back[panel]
        if self._ylim[panel] is None: self._ylim[panel] = b.get_ylim()
        elif self._ylim[panel] != b.get_ylim():
            raise RuntimeError("y limits changed after painting began; set them first")
        with np.errstate(all="ignore"):
            yp = b.transData.transform(np.column_stack([np.full(len(y), float(self.nmin)), y]))[:, 1]
        row = np.floor(yp - self.y0[panel])
        ok = np.isfinite(row) & (row >= 0) & (row < self.ph)
        return (n[ok] - self.nmin)*self.k, row[ok].astype(np.int64)

    def _over(self, panel, rr, cc, rgb, a):
        """Composite (rgb, alpha a) over the canvas at pixels (rr, cc); a may be per-pixel."""
        cv = self.canvas[panel]; a = np.broadcast_to(np.asarray(a, np.float32), rr.shape)[:, None]
        cv[rr, cc] = np.concatenate([rgb*a, a], axis=1) + cv[rr, cc]*(1 - a)

    def points(self, n, y, color, panel=0, h=3, alpha=1.0):
        """One mark per (n, y): exactly k px wide (n's own columns) and h px tall."""
        c, r = self._cols_rows(n, y, panel)
        dr = np.arange(h) - h//2
        rr = (r[:, None, None] + dr[None, :, None]).ravel()
        cc = np.broadcast_to(c[:, None, None] + np.arange(self.k)[None, None, :],
                             (len(c), h, self.k)).ravel()
        ok = (rr >= 0) & (rr < self.ph)
        self._over(panel, rr[ok], cc[ok], _hex_rgb(color), alpha)

    def density(self, n, y, lo="#dde1e6", hi="#3d434b", panel=0):
        """Every (n, y) sample binned into its pixel; colour runs lo -> hi with log(count).
        Returns the largest per-pixel count, for the legend."""
        c, r = self._cols_rows(n, y, panel)
        cnt = np.bincount(r*self.N + (c//self.k), minlength=self.ph*self.N).reshape(self.ph, self.N)
        rr, nn = np.nonzero(cnt)
        t = (np.log(cnt[rr, nn])/np.log(max(cnt.max(), 2)))[:, None].astype(np.float32)
        rgb = _hex_rgb(lo)*(1 - t) + _hex_rgb(hi)*t
        px = np.concatenate([rgb, np.ones((len(rr), 1), np.float32)], 1)   # opaque
        for d in range(self.k):
            self.canvas[panel][rr, nn*self.k + d] = px
        return int(cnt.max())

    def save(self, path):
        for p, b in enumerate(self.back):
            if b.get_xlim() != (self.nmin - 0.5, self.nmax + 0.5):
                raise RuntimeError("x limits were changed; the grid needs nmin-0.5..nmax+0.5")
            if self._ylim[p] is not None and b.get_ylim() != self._ylim[p]:
                raise RuntimeError("y limits changed after painting")
            cv = self.canvas[p]; a = cv[..., 3:4]
            rgba = np.concatenate([np.where(a > 0, cv[..., :3]/np.maximum(a, 1e-12), 0), a], 2)
            self.fig.figimage(np.clip(rgba, 0, 1), xo=self.x0, yo=self.y0[p], origin="lower",
                              zorder=1)
        self.fig.savefig(path, dpi=DPI)              # deliberately NOT bbox_inches="tight"
        if self.fig.canvas.get_width_height() != (self.W, self.H):
            raise RuntimeError(f"figure is {self.fig.canvas.get_width_height()}, "
                               f"expected {(self.W, self.H)}")
        return path
