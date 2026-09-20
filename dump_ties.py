"""
dump_ties.py -- export every tie point (p*>1/2) for given n to Parquet, for plotting.

    python dump_ties.py --n 100 [--n 200 ...] [--data data/] [--cusps cusps/cusps_all.csv]
    python dump_ties.py --n 100 --force          # rebuild even if the partition exists

Writes, per n:
    <data>/ties/n=NNNNN/part.parquet     every tie point 0<i<j<n with i+j>n
    <data>/cusps/n=NNNNN/part.parquet    the cusp subset, plus cusp-to-cusp distance columns
    <data>/manifest.csv                  one line per built partition

Stored columns (float64 unless noted; see cusps_data.py for everything derived from them):
    i, j                int16     n is the partition key
    pstar               tie point p*(n,i,j)
    ln_fi               log of the common mass f(i)=f(j) at p*, after normalisation
    E_minus_Ehalf       E(n,p*) - E(n,1/2)
    S_minus             left slope numerator; left slope = S_minus/(p* q*)
                        S_plus is NOT stored: it is S_minus + (j-i)f(i), and for most tie points
                        the kink (j-i)f(i) is many orders below S_minus, so a stored double S_plus
                        loses it entirely to cancellation (72 of 2401 rows at n=100 had
                        S_plus == S_minus bit-for-bit).  ln_fi determines it exactly instead.
    F3                  (n+i-j)(i+j-2np*) + (j-np*)
    is_cusp             bool      taken from the certified table, not from this script
    certified_by        str       'double'/'iv50'/... for cusps, '' otherwise
    gap_prev, gap_next  distance in p* to the adjacent tie points (NaN at the ends)
    rank_in_n           int32     position in p*-sorted order within this n
  cusps only, additionally:
    cusp_gap_prev/next, cusp_intervening_prev/next, nb_i, nb_j (nearest cusp by |p*| difference)

Precision: the masses are normalised by their own sum before E is accumulated.  Without that they
carry a shared relative error ~6e-13 from exp/lgamma, which leaves only ~3.6 digits on
E - E(1/2) at n=2000; normalised, ~6.3 digits.  E(n,1/2) is computed the same way.
is_cusp is never decided here -- the certified decision lives in cusps_fast.py.
"""
import argparse, csv, math, os, time
import numpy as np
from math import lgamma
from numba import njit

TINY = 1e-290
SCHEMA_VERSION = 1

@njit(cache=True)
def _dump_kernel(n, lnC, out_i, out_j, out_p, out_lnf, out_E, out_Sm, out_F3):
    """One row per tie point (i<j, i+j>n).  Returns the number written."""
    f = np.empty(n+1); w = np.empty(n+1, np.int64); c = 0
    for i in range(1, n):
        for j in range(max(i+1, n-i+1), n):
            m = j - i
            lnrho = (lnC[i] - lnC[j]) / m
            p = 1.0/(1.0 + np.exp(-lnrho)); q = 1.0 - p; rho = p/q
            md = int(np.floor((n+1)*p))
            if md > n: md = n
            f[md] = np.exp(lnC[md] + md*np.log(p) + (n-md)*np.log(q))
            for k in range(md, 0, -1):
                f[k-1] = f[k] * k / ((n-k+1.0)*rho)
            for k in range(md, n):
                f[k+1] = f[k] * rho*(n-k) / (k+1.0)
            for k in range(n+1):
                if f[k] < TINY: f[k] = 0.0
            f[j] = f[i]
            # normalise (Neumaier sum) so the shared exp/lgamma scale error cancels
            s = 0.0; comp = 0.0
            for k in range(n+1):
                t = s + f[k]
                if abs(s) >= abs(f[k]): comp += (s - t) + f[k]
                else:                   comp += (f[k] - t) + s
                s = t
            s = s + comp
            if s > 0.0:
                inv = 1.0/s
                for k in range(n+1): f[k] *= inv
            sp_ = md
            if sp_ >= j: sp_ = j - 1
            if sp_ < i: sp_ = i
            a = 0; b = n; r = 0
            while a <= sp_ or b > sp_:
                if a > sp_: take_left = False
                elif b <= sp_: take_left = True
                elif f[a] < f[b]: take_left = True
                elif f[a] > f[b]: take_left = False
                else: take_left = not (a == i and b == j)
                if take_left: k = a; a += 1
                else:         k = b; b -= 1
                w[k] = r; r += 1
            Sm = 0.0; E = 0.0; ce = 0.0
            for k in range(n+1):
                t = w[k]*f[k]
                Sm += t*(k - n*p)
                u = E + t                      # Neumaier again for E
                if abs(E) >= abs(t): ce += (E - u) + t
                else:                ce += (t - u) + E
                E = u
            E = E + ce
            out_i[c] = i; out_j[c] = j; out_p[c] = p
            out_lnf[c] = lnC[i] + i*np.log(p) + (n-i)*np.log(q) - np.log(s)
            out_E[c] = E; out_Sm[c] = Sm
            out_F3[c] = (n+i-j)*(i+j-2*n*p) + (j-n*p)
            c += 1
    return c

def lnC_arr(n):
    return np.array([lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1) for k in range(n+1)])

def E_half(n):
    """E(n,1/2) with the same normalisation convention."""
    lnC = lnC_arr(n)
    f = np.exp(lnC - n*math.log(2.0))
    f = f/math.fsum(f.tolist())
    g = np.sort(f)
    return math.fsum((np.arange(n+1)*g).tolist())

def n_ties(n):
    return sum(max(0, (n-1) - max(i+1, n-i+1) + 1) for i in range(1, n))

def certified_cusps(path, n):
    """{(i,j): certified_by} for this n, from the certified table."""
    out = {}
    if not os.path.exists(path): return out
    with open(path) as fh:
        for r in csv.DictReader(fh):
            if int(r['n']) == n and r['pstar'] != '':
                out[(int(r['i']), int(r['j']))] = r['certified_by']
    return out

def build(n, data, cusp_csv, force=False):
    import pyarrow as pa, pyarrow.parquet as pq
    tdir = os.path.join(data, "ties", f"n={n:05d}"); cdir = os.path.join(data, "cusps", f"n={n:05d}")
    if os.path.exists(os.path.join(tdir, "part.parquet")) and not force:
        print(f"n={n}: exists, skipping (use --force)"); return None
    t0 = time.time()
    cap = n_ties(n); lnC = lnC_arr(n)
    oi = np.empty(cap, np.int64); oj = np.empty(cap, np.int64)
    op = np.empty(cap); olf = np.empty(cap); oE = np.empty(cap)
    oSm = np.empty(cap); oF3 = np.empty(cap)
    c = _dump_kernel(n, lnC, oi, oj, op, olf, oE, oSm, oF3)
    assert c == cap, f"kernel wrote {c} rows, expected {cap}"
    Eh = E_half(n)
    o = np.argsort(op, kind='stable')
    oi, oj, op, olf, oE, oSm, oF3 = (x[o] for x in (oi, oj, op, olf, oE, oSm, oF3))
    gap_next = np.full(c, np.nan); gap_prev = np.full(c, np.nan)
    if c > 1:
        d = np.diff(op); gap_next[:-1] = d; gap_prev[1:] = d
    cusps = certified_cusps(cusp_csv, n)
    is_cusp = np.array([(int(a), int(b)) in cusps for a, b in zip(oi, oj)])
    cert = np.array([cusps.get((int(a), int(b)), '') for a, b in zip(oi, oj)])
    cols = dict(i=pa.array(oi.astype(np.int16)), j=pa.array(oj.astype(np.int16)),
                pstar=pa.array(op), ln_fi=pa.array(olf), E_minus_Ehalf=pa.array(oE - Eh),
                S_minus=pa.array(oSm), F3=pa.array(oF3),
                is_cusp=pa.array(is_cusp), certified_by=pa.array(cert).dictionary_encode(),
                gap_prev=pa.array(gap_prev), gap_next=pa.array(gap_next),
                rank_in_n=pa.array(np.arange(c, dtype=np.int32)))
    meta = {b'n': str(n).encode(), b'E_half': repr(Eh).encode(), b'n_ties': str(c).encode(),
            b'n_cusps': str(int(is_cusp.sum())).encode(), b'normalised_masses': b'1',
            b'schema_version': str(SCHEMA_VERSION).encode(), b'TINY': repr(TINY).encode()}
    def write(table, d):
        os.makedirs(d, exist_ok=True)
        table = table.replace_schema_metadata(meta)
        pq.write_table(table, os.path.join(d, "part.parquet"), compression='zstd',
                       use_dictionary=False,
                       use_byte_stream_split=[k for k in table.column_names
                                              if pa.types.is_floating(table.schema.field(k).type)])
    write(pa.table(cols), tdir)
    # --- cusp subset with cusp-to-cusp distances ---
    ci = np.flatnonzero(is_cusp)
    sub = {k: v.take(pa.array(ci)) for k, v in cols.items()}
    cp = op[ci]
    cg_prev = np.full(len(ci), np.nan); cg_next = np.full(len(ci), np.nan)
    iv_prev = np.full(len(ci), -1, np.int32); iv_next = np.full(len(ci), -1, np.int32)
    if len(ci) > 1:
        d = np.diff(cp); cg_next[:-1] = d; cg_prev[1:] = d
        gapidx = np.diff(ci) - 1
        iv_next[:-1] = gapidx; iv_prev[1:] = gapidx
    near = np.where(np.nan_to_num(cg_prev, nan=np.inf) <= np.nan_to_num(cg_next, nan=np.inf),
                    np.arange(len(ci))-1, np.arange(len(ci))+1)
    near = np.clip(near, 0, max(len(ci)-1, 0))
    sub.update(cusp_gap_prev=pa.array(cg_prev), cusp_gap_next=pa.array(cg_next),
               cusp_intervening_prev=pa.array(iv_prev), cusp_intervening_next=pa.array(iv_next),
               nb_i=pa.array(oi[ci][near].astype(np.int16)), nb_j=pa.array(oj[ci][near].astype(np.int16)))
    write(pa.table(sub), cdir)
    dt = time.time()-t0
    line = dict(n=n, n_ties=c, n_cusps=int(is_cusp.sum()), E_half=repr(Eh),
                schema_version=SCHEMA_VERSION, seconds=round(dt, 2),
                built=time.strftime("%Y-%m-%dT%H:%M:%S"))
    mf = os.path.join(data, "manifest.csv"); new = not os.path.exists(mf)
    with open(mf, 'a', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(line))
        if new: w.writeheader()
        w.writerow(line)
    print(f"n={n}: {c:,} ties, {int(is_cusp.sum())} cusps, {dt:.1f}s -> {tdir}")
    return line

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, action='append', required=True)
    ap.add_argument("--data", default="data")
    ap.add_argument("--cusps", default="cusps/cusps_all.csv")
    ap.add_argument("--force", action='store_true')
    a = ap.parse_args()
    for n in a.n: build(n, a.data, a.cusps, a.force)
