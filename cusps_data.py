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

# ---------------------------------------------------------------- reading the shared store
#
# Since 2026-10-12 everything here comes from OBD-core's shared store (obd_core.store, ~/obd-data,
# every value from obd_core.tie_table).  The ``data`` arguments are kept so old calls still work, and
# are ignored.  The columns are built exactly as dump_ties.py built its Parquet files, which these
# functions read before; old and new were compared column by column, bit for bit.

_STORE = None


def _store():
    global _STORE
    if _STORE is None:
        from obd_core import store
        _STORE = store.Store()
    return _STORE


def _run_end(kind):
    """The last n of the store's unbroken run of complete tables of this kind, from n = 3."""
    st, n = _store(), 2
    while st.covered(kind, n + 1) == [[0.5, 1.0]]:
        n += 1
    return n


def _dump_ties_columns(n, r):
    """dump_ties.py's columns for one n, from tie_table(n) rows (p* >= 1/2, the axis as row 0):
    {"ties": every tie point, "cusps": the cusp rows with their cusp-to-cusp extras}."""
    import obd_core as core
    is_cusp, decided, n_pairs = r["is_cusp"], np.asarray(r["decided_by"], dtype=object), r["n_tied_pairs"]
    c = len(r["i"])
    gap_next = np.full(c, np.nan); gap_prev = np.full(c, np.nan)
    if c > 1:
        d = np.diff(r["pstar"]); gap_next[:-1] = d; gap_prev[1:] = d
    cols = dict(i=r["i"].astype(np.int16), j=r["j"].astype(np.int16), pstar=r["pstar"], ln_fi=r["ln_fi"],
                E_minus_Ehalf=r["E"] - core.E_half(n), S_minus=r["S_minus"], F3=r["F3"], is_cusp=is_cusp,
                decided_by=decided, gap_prev=gap_prev, gap_next=gap_next,
                rank_in_n=np.arange(c, dtype=np.int32), n_tied_pairs=n_pairs.astype(np.int16))
    ci = np.flatnonzero(is_cusp)
    sub = {k: v[ci] for k, v in cols.items()}
    cp = r["pstar"][ci]
    cg_prev = np.full(len(ci), np.nan); cg_next = np.full(len(ci), np.nan)
    iv_prev = np.full(len(ci), -1, np.int32); iv_next = np.full(len(ci), -1, np.int32)
    if len(ci) > 1:
        d = np.diff(cp); cg_next[:-1] = d; cg_prev[1:] = d
        gi = np.diff(ci) - 1; iv_next[:-1] = gi; iv_prev[1:] = gi
    near = np.where(np.nan_to_num(cg_prev, nan=np.inf) <= np.nan_to_num(cg_next, nan=np.inf),
                    np.arange(len(ci)) - 1, np.arange(len(ci)) + 1)
    near = np.clip(near, 0, max(len(ci) - 1, 0))
    sub.update(cusp_gap_prev=cg_prev, cusp_gap_next=cg_next, cusp_intervening_prev=iv_prev,
               cusp_intervening_next=iv_next, nb_i=r["i"][ci][near].astype(np.int16),
               nb_j=r["j"][ci][near].astype(np.int16))
    return {"ties": cols, "cusps": sub}


def _read(kind, n=None, data="data", columns=None):
    import obd_core as core
    avail = available(kind)
    if n is None:
        ns = avail
    elif isinstance(n, tuple):
        ns = [v for v in avail if n[0] <= v <= n[1]]
    else:
        ns = [int(n)]                       # any n: the store reads it, or computes and saves it
    if not ns:
        raise FileNotFoundError(f"no n matches {n} in the store's tie tables")
    out = {}; nn = []; eh = []
    for v in ns:
        t = _dump_ties_columns(v, _store().ties(v))[kind]
        for k in (columns or t):
            out.setdefault(k, []).append(t[k])
        rows = len(t["pstar"])
        nn.append(np.full(rows, v, np.int32)); eh.append(np.full(rows, core.E_half(v)))
    d = {k: np.concatenate(v) for k, v in out.items()}
    d["n"] = np.concatenate(nn)
    d["E_half"] = np.concatenate(eh)
    return d


CATALOGUE = ["n", "i", "j", "pstar", "E", "F3", "F3_sign", "S_minus", "S_plus", "slope_left",
             "slope_right", "certified_by"]


def read_cusps(columns=None, nmin=3, nmax=None, exact_kink=False, as_printed=True):
    """The cusp catalogue that cusps/cusps_all.csv held, as a pyarrow Table, from the shared store.

    One row per certified cusp with p* > 1/2 (the axis row excluded), for n = nmin .. nmax (default:
    to the end of the store's unbroken run of complete cusp tables), in (n, i, j) order as the CSV
    had them, with the CSV's columns and types (n, i, j int64; F3_sign "+"/"-"; certified_by = the
    store's decided_by).  pstar, E, F3, S_minus and slope_left are bit for bit the values the CSV
    was printed from (obd_core.evaluate).  S_plus = S_minus + (j-i) f(i) and slope_right =
    S_plus/(p q), as the CSV defined them; f(i) comes from the stored ln f(i), so these two can
    differ from the CSV's in the last bit or two -- exact_kink=True takes them from
    obd_core.evaluate instead (exact; about 30 us a row, 2-3 minutes for the whole catalogue).  The cusp sets are identical.
    as_printed (default) rounds every float through "%.16g", as the CSV printed it, so the values
    are bit for bit what reading the CSV gave (with exact_kink, for every column).
    Replaces pc.read_csv("cusps/cusps_all.csv", ...): pass the same column list."""
    import pyarrow as pa
    import obd_core as core
    st = _store()
    cols = list(columns or CATALOGUE)
    need = {"pstar", "i", "j", "F3", "S_minus", "ln_fi", "decided_by", "E"}
    hi = _run_end("cusps") if nmax is None else int(nmax)
    acc = {c: [] for c in cols}
    for n in range(int(nmin), hi + 1):
        t = st.cusps(n, columns=sorted(need))
        k = np.flatnonzero(t["pstar"] > 0.5)
        k = k[np.lexsort((t["j"][k], t["i"][k]))]
        g = {name: v[k] for name, v in t.items()}
        pq_ = g["pstar"] * (1 - g["pstar"])
        if exact_kink:
            lnC = core.lnC_arr(n)
            Sp = np.array([core.evaluate(n, int(a), int(b), lnC)[4] for a, b in zip(g["i"], g["j"])], float)
        else:
            Sp = g["S_minus"] + (g["j"] - g["i"]) * np.exp(g["ln_fi"])
        derived = {"n": np.full(k.size, n, np.int64), "i": g["i"].astype(np.int64), "j": g["j"].astype(np.int64),
                   "F3_sign": np.where(g["F3"] > 0, "+", "-").astype(object), "S_plus": Sp,
                   "slope_left": g["S_minus"] / pq_, "slope_right": Sp / pq_,
                   "certified_by": np.asarray(g["decided_by"], dtype=object)}
        for c in cols:
            acc[c].append(derived[c] if c in derived else g[c])
    out = {}
    for c in cols:
        v = np.concatenate(acc[c]) if acc[c] else np.empty(0)
        if as_printed and v.dtype.kind == "f":
            v = np.char.mod("%.16g", v).astype(np.float64)
        out[c] = pa.array(v)
    return pa.table(out)


def cusp_lines(nmin=3, nmax=None, exact_kink=False):
    """The catalogue as the text lines cusps_fast.py wrote (header first), each value printed with
    the same "%.16g" -- a drop-in for reading cusps/cusps_all.csv as text."""
    t = read_cusps(None, nmin, nmax, exact_kink, as_printed=False)
    yield ",".join(CATALOGUE) + "\n"
    cols = [t[c].to_numpy(zero_copy_only=False) for c in CATALOGUE]
    for r in zip(*cols):
        n, i, j, p, E, F3, sg, Sm, Sp, sl, sr, how = r
        yield (f"{n},{i},{j},{p:.16g},{E:.16g},{F3:.16g},{sg},{Sm:.16g},{Sp:.16g},{sl:.16g},{sr:.16g},{how}\n")


def cusp_rows(nmin=3, nmax=None, exact_kink=False):
    """The catalogue as csv.DictReader rows (strings, formatted as the CSV had them)."""
    import csv
    return csv.DictReader(cusp_lines(nmin, nmax, exact_kink))


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
    """The n whose full tie tables the store holds (both kinds here are built from those)."""
    st = _store()
    return [n for n in st.ns("ties") if st.covered("ties", n) == [[0.5, 1.0]]]


def catalogue(path=None, columns=None, exact_kink=False, nmax=None):
    """read_cusps(columns), or, given a path, that CSV file read as before (pyarrow)."""
    if path is None:
        return read_cusps(columns, nmax=nmax, exact_kink=exact_kink)
    import pyarrow.csv as pc
    opts = pc.ConvertOptions(include_columns=list(columns)) if columns else None
    return pc.read_csv(path, convert_options=opts)


def catalogue_rows(path=None, exact_kink=False, nmax=None):
    """cusp_rows(), or, given a path, that CSV file through csv.DictReader as before."""
    if path is None:
        return cusp_rows(nmax=nmax, exact_kink=exact_kink)
    import csv
    return csv.DictReader(open(path))
