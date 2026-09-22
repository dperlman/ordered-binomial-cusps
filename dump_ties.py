"""
dump_ties.py -- export every tie point (p*>1/2) for given n to Parquet, for plotting.

    python dump_ties.py --n 100 [--n 200 ...] [--data data/] [--force] [--workers 8]
    python dump_ties.py --n 100 --verify cusps/cusps_all.csv     # cross-check against the CSV table

Writes, per n:
    <data>/ties/n=NNNNN/part.parquet     every tie point 0<=i<j<=n with i+j>n
    <data>/cusps/n=NNNNN/part.parquet    the cusp subset, plus cusp-to-cusp distance columns
    <data>/manifest.csv                  one line per built partition

The mathematics is binom_core.py -- the same kernel and the same certify() the certified generator
uses.  Certification is redone HERE rather than looked up: a lookup silently mislabels every tie
point as a non-cusp for any n the CSV table does not cover (e.g. n=3000 today), which looks like a
result instead of an error.  Re-certifying costs only the CHECK-tagged tie points, 4-86 per n
(~17 s at n=2000 against ~23 s of screening).  --verify compares the outcome against the certified
CSV where it exists, which is then a regression test rather than a data dependency.

Stored columns (float64 unless noted; see cusps_data.py for everything derived from them):
    i, j                int16     n is the partition key
    pstar               tie point p*(n,i,j)
    ln_fi               log of the common mass f(i)=f(j) at p*, after normalisation
    E_minus_Ehalf       E(n,p*) - E(n,1/2)
    S_minus             left slope numerator; left slope = S_minus/(p* q*).
                        S_plus is NOT stored: it is S_minus + (j-i)f(i), and for most tie points
                        the kink is many orders below S_minus, so a stored double S_plus loses it
                        to cancellation (72 of 2401 rows at n=100 had S_plus == S_minus exactly).
                        ln_fi determines it exactly instead.
    F3                  (n+i-j)(i+j-2np*) + (j-np*)
    is_cusp             bool
    decided_by          str       'double', or 'iv50'/'iv100'/'iv200' when interval arithmetic was
                                  needed; 'UNRESOLVED' if even 200 digits did not settle it
    gap_prev, gap_next  distance in p* to the adjacent tie points (NaN at the ends)
    rank_in_n           int32     position in p*-sorted order within this n
    n_tied_pairs        int16     1 for an ordinary tie point.  Row 0 of every partition is the
                                  SYMMETRY AXIS p=1/2, where all mirror pairs (i,n-i) tie at once;
                                  there it counts them, and (i,j) is the sentinel (0,n) so that
                                  band = i+j-n = 0 marks it (every real tie point has band >= 1).
                                  E_minus_Ehalf is exactly 0 and u = S_-/kappa is exactly -1/2
                                  there, by E(p)=E(1-p).  F3 is NaN: the pair decomposition does
                                  not apply to a multi-tie.
  cusps only, additionally:
    cusp_gap_prev/next, cusp_intervening_prev/next, nb_i, nb_j (nearest cusp by |p*| difference)

Masses are normalised by their own sum before E is accumulated; see binom_core.py for why.
"""
import argparse, csv, os, time
from multiprocessing import Pool
import numpy as np

import binom_core as core
from binom_core import TAG_CHECK, TAG_MIN

SCHEMA_VERSION = 3

def _screen_chunk(a):
    n, lo, hi = a
    return core.screen(n, collect_all=True, i_lo=lo, i_hi=hi)

def _certify_one(a):
    n, i, j = a
    v, how = core.certify_escalating(n, i, j)
    return (v == 'MIN'), (how if v else 'UNRESOLVED')

def build(n, data, force=False, verify=None, workers=1, pool=None):
    import pyarrow as pa, pyarrow.parquet as pq
    tdir = os.path.join(data, "ties", f"n={n:05d}"); cdir = os.path.join(data, "cusps", f"n={n:05d}")
    if os.path.exists(os.path.join(tdir, "part.parquet")) and not force:
        print(f"n={n}: exists, skipping (use --force)"); return None
    t0 = time.time()
    r = _screen_parallel(n, workers, pool)
    for k in ('tag_alt', 'rbnd'): r.pop(k, None)   # trigger diagnostics; not part of the dataset
    o = np.argsort(r['pstar'], kind='stable')
    r = {k: v[o] for k, v in r.items()}
    # prepend the symmetry axis p=1/2 (see binom_core.axis_point)
    aSm, aSp, aE, akap, apairs = core.axis_point(n)
    axis = dict(i=0, j=n, pstar=0.5, ln_fi=np.log(akap/n), E=aE, S_minus=aSm,
                F3=np.nan, tag=(TAG_MIN if aSm < 0 < aSp else 0))
    for k in r: r[k] = np.concatenate([np.array([axis[k]], dtype=r[k].dtype), r[k]])
    n_pairs = np.ones(len(r['i']), np.int16); n_pairs[0] = apairs
    c = len(r['i'])
    # --- certify every CHECK-tagged tie point, exactly as the generator does ---
    decided = np.array(['double']*c, dtype=object)
    decided[0] = 'symmetry'            # the axis row is settled exactly by E(p) = E(1-p)
    is_cusp = r['tag'] == TAG_MIN
    checks = np.flatnonzero(r['tag'] == TAG_CHECK)
    t_cert = time.time()
    if len(checks):
        args = [(n, int(r['i'][t]), int(r['j'][t])) for t in checks]
        res = pool.map(_certify_one, args, chunksize=1) if pool else [_certify_one(a) for a in args]
        for t, (cusp, how) in zip(checks, res):
            is_cusp[t] = cusp; decided[t] = how
    t_cert = time.time() - t_cert
    n_unres = int((decided == 'UNRESOLVED').sum())
    Eh = core.E_half(n)
    gap_next = np.full(c, np.nan); gap_prev = np.full(c, np.nan)
    if c > 1:
        d = np.diff(r['pstar']); gap_next[:-1] = d; gap_prev[1:] = d
    cols = dict(i=pa.array(r['i'].astype(np.int16)), j=pa.array(r['j'].astype(np.int16)),
                pstar=pa.array(r['pstar']), ln_fi=pa.array(r['ln_fi']),
                E_minus_Ehalf=pa.array(r['E'] - Eh), S_minus=pa.array(r['S_minus']),
                F3=pa.array(r['F3']), is_cusp=pa.array(is_cusp),
                decided_by=pa.array(decided.tolist()).dictionary_encode(),
                gap_prev=pa.array(gap_prev), gap_next=pa.array(gap_next),
                rank_in_n=pa.array(np.arange(c, dtype=np.int32)),
                n_tied_pairs=pa.array(n_pairs))
    meta = {b'n': str(n).encode(), b'E_half': repr(Eh).encode(), b'n_ties': str(c).encode(),
            b'n_cusps': str(int(is_cusp.sum())).encode(), b'n_checked': str(len(checks)).encode(),
            b'n_unresolved': str(n_unres).encode(), b'normalised_masses': b'1',
            b'includes_axis_p_half': b'1', b'axis_n_tied_pairs': str(apairs).encode(),
            b'schema_version': str(SCHEMA_VERSION).encode(), b'TINY': repr(core.TINY).encode(),
            b'MARGIN': repr(core.MARGIN).encode(), b'GAP': repr(core.GAP).encode()}
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
    cp = r['pstar'][ci]
    cg_prev = np.full(len(ci), np.nan); cg_next = np.full(len(ci), np.nan)
    iv_prev = np.full(len(ci), -1, np.int32); iv_next = np.full(len(ci), -1, np.int32)
    if len(ci) > 1:
        d = np.diff(cp); cg_next[:-1] = d; cg_prev[1:] = d
        gi = np.diff(ci) - 1; iv_next[:-1] = gi; iv_prev[1:] = gi
    near = np.where(np.nan_to_num(cg_prev, nan=np.inf) <= np.nan_to_num(cg_next, nan=np.inf),
                    np.arange(len(ci))-1, np.arange(len(ci))+1)
    near = np.clip(near, 0, max(len(ci)-1, 0))
    sub.update(cusp_gap_prev=pa.array(cg_prev), cusp_gap_next=pa.array(cg_next),
               cusp_intervening_prev=pa.array(iv_prev), cusp_intervening_next=pa.array(iv_next),
               nb_i=pa.array(r['i'][ci][near].astype(np.int16)),
               nb_j=pa.array(r['j'][ci][near].astype(np.int16)))
    write(pa.table(sub), cdir)
    dt = time.time()-t0
    msg = (f"n={n}: {c:,} rows (incl. p=1/2 axis, {apairs} pairs), {int(is_cusp.sum())} cusps, "
           f"{len(checks)} certified "
           f"({t_cert:.1f}s of {dt:.1f}s)")
    if n_unres: msg += f"  *** {n_unres} UNRESOLVED ***"
    print(msg + f" -> {tdir}")
    if verify: verify_against_csv(n, verify, r, is_cusp, decided)
    line = dict(n=n, n_ties=c, n_cusps=int(is_cusp.sum()), n_checked=len(checks),
                n_unresolved=n_unres, E_half=repr(Eh), schema_version=SCHEMA_VERSION,
                seconds=round(dt, 2), built=time.strftime("%Y-%m-%dT%H:%M:%S"))
    mf = os.path.join(data, "manifest.csv"); new = not os.path.exists(mf)
    with open(mf, 'a', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(line))
        if new: w.writeheader()
        w.writerow(line)
    return line

def _screen_parallel(n, workers, pool):
    """Screen one n, splitting the i-loop into work-balanced chunks across `pool`.

    Each tie point is computed identically however the range is cut, and the chunks are
    concatenated in i order, so the result is independent of `workers`.
    """
    if not pool or workers <= 1:
        return core.screen(n, collect_all=True)
    chunks = core.work_chunks(n, workers)
    parts = pool.map(_screen_chunk, [(n, a, b) for a, b in chunks])
    return {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}

def verify_against_csv(n, path, r, is_cusp, decided):
    """Regression check: do we agree with the certified CSV table where it covers this n?"""
    if not os.path.exists(path):
        print(f"  verify: {path} not found, skipped"); return
    ref = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            if int(row['n']) == n and row['pstar'] != '':
                ref[(int(row['i']), int(row['j']))] = row['certified_by']
    if not ref:
        print(f"  verify: the CSV table has no rows for n={n} "
              f"(this is exactly the case a lookup would have silently mislabelled)"); return
    got = {(int(a), int(b)) for a, b, c in zip(r['i'], r['j'], is_cusp) if c}
    got.discard((0, n))                      # the axis row is not in the certified CSV table
    missing = set(ref) - got; extra = got - set(ref)
    same_route = sum(1 for t in range(len(r['i']))
                     if (int(r['i'][t]), int(r['j'][t])) in ref
                     and ref[(int(r['i'][t]), int(r['j'][t]))] == decided[t])
    if missing or extra:
        print(f"  verify: *** MISMATCH *** missing {sorted(missing)[:5]} extra {sorted(extra)[:5]}")
    else:
        print(f"  verify: cusp set identical to the certified table ({len(ref)} cusps); "
              f"certification route agrees on {same_route}/{len(ref)}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, action='append', required=True)
    ap.add_argument("--data", default="data")
    ap.add_argument("--force", action='store_true')
    ap.add_argument("--workers", type=int, default=8,
                    help="processes for the i-loop and the certifications (1 = serial)")
    ap.add_argument("--verify", nargs='?', const="cusps/cusps_all.csv", default=None,
                    help="cross-check the cusp set against this certified CSV table")
    a = ap.parse_args()
    core.screen(10, collect_all=True)     # compile the kernel once up front
    if a.workers > 1:
        with Pool(a.workers) as pool:
            for n in a.n: build(n, a.data, a.force, a.verify, a.workers, pool)
    else:
        for n in a.n: build(n, a.data, a.force, a.verify)
