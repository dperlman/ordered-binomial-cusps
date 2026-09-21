"""
cusps_data.py -- read the Parquet tie-point / cusp datasets and derive the plotting quantities.

    import cusps_data as cd
    df = cd.load_ties(100)                 # one n, all tie points   (pandas-free: dict of arrays)
    df = cd.load_cusps()                   # all n present, cusp rows only
    df = cd.load_cusps(n=(100, 2000))      # a range of n
    cd.columns()                           # what is stored vs derived

Only pstar, ln_fi, E_minus_Ehalf, S_minus, S_plus and F3 are stored; everything below is derived
here so the two can never drift apart.  Conventions follow RESEARCH_LOG.md:

    E = B + V,  B = A + T          B smooth backbone, V = (1/2)|f_j - f_i| the kink (V=0 at p*)
    w_i = n+i-j = r,  f = f(i) = f(j) at p*,  slope_unit = f/(p* q*)
    kink_pos = S_-/kappa   -- where zero falls inside the slope jump;  cusp <=> -1 < kink_pos < 0,
               and kink_pos = -1/2 exactly on the p=1/2 axis row, by E(p) = E(1-p)
    T(p*)  = (2r+1) f          V(p*) = 0          A(p*) = E - T   (all NaN on the axis row)
    slopes:  E'_- = S_minus/(p q)         E'_+ = S_plus/(p q)
             (T+V)'_+ = F3 * slope_unit              (T+V)'_- = (F3-(j-i)) * slope_unit
             T'       = (F3-(j-i)/2) * slope_unit    V'_±     = +-(j-i)/2 * slope_unit
             A'       = E'_+ - F3 * slope_unit       (smooth; equals E'_- - (F3-(j-i)) u)
             D        = (j-i) * slope_unit           the slope jump, = E'_+ - E'_-
S_plus is derived as S_minus + (j-i)f, and slope_right as slope_left + D.  Do not compute
slope_right as S_plus/(p q): for most tie points the kink is orders of magnitude below S_minus and
the subtraction loses it completely.
F3 is therefore the right-hand slope of T+V in units of slope_unit, NOT the slope of T alone.
"""
import glob, os
import numpy as np

STORED = ["i", "j", "pstar", "ln_fi", "E_minus_Ehalf", "S_minus", "F3",
          "is_cusp", "decided_by", "gap_prev", "gap_next", "rank_in_n", "n_tied_pairs"]
CUSP_EXTRA = ["cusp_gap_prev", "cusp_gap_next", "cusp_intervening_prev", "cusp_intervening_next",
              "nb_i", "nb_j"]
DERIVED = ["n", "is_axis", "width", "band", "w_i", "f_i", "slope_unit", "kink_pos", "kappa",
           "S_plus", "E", "E_half", "T", "A", "V",
           "slope_left", "slope_right", "slope_T", "slope_A", "slope_V_right", "D",
           "gap_nearest", "n_gap_nearest", "cusp_gap_nearest", "is_first_band"]

def columns():
    return {"stored": STORED, "stored (cusp files only)": CUSP_EXTRA, "derived here": DERIVED}

def _read(kind, n=None, data="data", columns=None):
    import pyarrow.parquet as pq, pyarrow.dataset as ds
    root = os.path.join(data, kind)
    parts = sorted(glob.glob(os.path.join(root, "n=*")))
    if not parts: raise FileNotFoundError(f"no partitions under {root} (build with dump_ties.py)")
    def keep(p):
        v = int(os.path.basename(p).split("=")[1])
        if n is None: return True
        if isinstance(n, tuple): return n[0] <= v <= n[1]
        return v == n
    parts = [p for p in parts if keep(p)]
    if not parts: raise FileNotFoundError(f"no partition matches n={n} under {root}")
    out = {}; ns = []; ehalf = {}
    for p in parts:
        fn = os.path.join(p, "part.parquet")
        pf = pq.ParquetFile(fn)
        md = pf.schema_arrow.metadata or {}
        nn = int(os.path.basename(p).split("=")[1])
        ehalf[nn] = float(md.get(b"E_half", b"nan"))
        t = pf.read(columns=columns)
        for k in t.column_names:
            out.setdefault(k, []).append(t.column(k).to_numpy(zero_copy_only=False))
        ns.append(np.full(t.num_rows, nn, np.int32))
    d = {k: np.concatenate(v) for k, v in out.items()}
    d["n"] = np.concatenate(ns)
    d["E_half"] = np.array([ehalf[x] for x in d["n"]])
    return d

def _derive(d):
    n = d["n"].astype(np.float64); i = d["i"].astype(np.float64); j = d["j"].astype(np.float64)
    p = d["pstar"]; q = 1.0 - p; pq_ = p*q
    d["is_axis"] = d["n_tied_pairs"] > 1 if "n_tied_pairs" in d else np.zeros(len(p), bool)
    d["width"] = (j - i).astype(np.int32)
    d["band"] = (i + j - n).astype(np.int32)
    d["is_first_band"] = d["band"] == 1
    d["w_i"] = (n + i - j).astype(np.int32)
    d["f_i"] = np.exp(d["ln_fi"])
    d["slope_unit"] = d["f_i"]/pq_          # the unit f/(p*q*) that the slopes below are in
    d["E"] = d["E_half"] + d["E_minus_Ehalf"]
    r = d["w_i"].astype(np.float64)
    d["T"] = (2*r + 1)*d["f_i"]
    d["V"] = np.zeros_like(p)
    d["A"] = d["E"] - d["T"]
    d["kappa"] = (j - i)*d["f_i"]               # the kink (j-i)f(i); S_plus - S_minus
    d["S_plus"] = d["S_minus"] + d["kappa"]
    with np.errstate(divide="ignore", invalid="ignore"):
        d["kink_pos"] = d["S_minus"]/d["kappa"]   # where zero sits in the jump; cusp <=> -1<.<0
    d["slope_left"] = d["S_minus"]/pq_
    d["D"] = d["kappa"]/pq_
    d["slope_right"] = d["slope_left"] + d["D"]  # never S_plus/pq: that cancels the kink away
    d["slope_T"] = (d["F3"] - (j - i)/2.0)*d["slope_unit"]
    d["slope_V_right"] = (j - i)/2.0*d["slope_unit"]
    d["slope_A"] = d["slope_right"] - d["F3"]*d["slope_unit"]
    for a, b, out in (("gap_prev", "gap_next", "gap_nearest"),
                      ("cusp_gap_prev", "cusp_gap_next", "cusp_gap_nearest")):
        if a in d:
            d[out] = np.fmin(np.nan_to_num(d[a], nan=np.inf), np.nan_to_num(d[b], nan=np.inf))
            d[out][np.isinf(d[out])] = np.nan
    d["n_gap_nearest"] = d["n"]*d["gap_nearest"]
    # the pair decomposition T/A/V is meaningless on the p=1/2 axis row (it is a multi-tie), but
    # u, kappa, D and the slopes are all still correct there.
    if d["is_axis"].any():
        for c in ("T", "A", "V", "slope_T", "slope_V_right", "slope_A", "w_i", "f_i"):
            d[c] = np.where(d["is_axis"], np.nan, d[c])
    return d

def load_ties(n, data="data", columns=None):
    return _derive(_read("ties", n, data, columns))

def load_cusps(n=None, data="data", columns=None):
    return _derive(_read("cusps", n, data, columns))

def available(kind="ties", data="data"):
    return sorted(int(os.path.basename(p).split("=")[1])
                  for p in glob.glob(os.path.join(data, kind, "n=*")))
