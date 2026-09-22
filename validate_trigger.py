"""validate_trigger.py -- regression test for the sharpened CHECK trigger (binom_core).

Two directions, both needed:

  --logged  Every tie point the OLD rule sent to interval arithmetic, read from
            cusps/interval_checks.log (n,i,j,verdict,route).  For each one, confirm the sharpened
            trigger decides it in double precision AND that the double verdict equals the verdict
            mpmath certified.  This is the "is anything the new rule drops actually wrong?" test,
            and it is free: the certified answers are already on record.

  --sweep   Every tie point of the given n, comparing both tags.  This is the other direction --
            can the new bound flag something the old rule decided?  A near-mode cluster with
            f ~ 1e-2 and |a| ~ n would give a bound of order 10 and must raise the check count.

Run: .venv/bin/python validate_trigger.py --logged
     .venv/bin/python validate_trigger.py --sweep 1100 1500 2000 3000
"""
import argparse, collections, sys
import numpy as np
import binom_core as core

C, MIN = core.TAG_CHECK, core.TAG_MIN

def logged(path):
    rows = collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            n, i, j, verdict = line.split(',')[:4]
            rows[int(n)].append((int(i), int(j), verdict))
    tot = sum(len(v) for v in rows.values())
    print(f"{tot} logged interval checks over {len(rows)} values of n "
          f"({min(rows)}..{max(rows)})")
    still, wrong, ok = [], [], 0
    for n in sorted(rows):
        lnC = core.lnC_arr(n)
        f = np.empty(n+1); w = np.empty(n+1, np.int64)
        for i, j, verdict in rows[n]:
            *_, tag_old, tag_new, rbnd = core._one_tie(n, lnC, i, j, f, w)
            if tag_old != C:
                wrong.append((n, i, j, 'old rule did not flag this at all'))
            if tag_new == C:
                still.append((n, i, j, verdict, rbnd))
            elif ('MIN' if tag_new == MIN else 'NOT') != verdict:
                wrong.append((n, i, j, f"double says {'MIN' if tag_new==MIN else 'NOT'}, "
                                       f"certified {verdict}"))
            else:
                ok += 1
    print(f"  decided in double by the sharpened trigger, verdict MATCHES certified: {ok}")
    print(f"  still CHECK (would still go to mpmath):                                {len(still)}")
    print(f"  DISAGREEMENTS:                                                         {len(wrong)}")
    for r in wrong[:20]: print("    ", r)
    for r in still[:20]: print("    still: ", r)
    return len(wrong)

def sweep(ns):
    print(f"{'n':>6} {'ties':>10} {'old':>8} {'new':>8} {'dropped':>8} {'added':>7} "
          f"{'maxR dropped':>13} {'maxR all':>11}")
    bad = 0
    for n in ns:
        r = core.screen(n, collect_all=True, sharp=False)   # so tag=old, tag_alt=new
        told, tnew, rb = r['tag'], r['tag_alt'], r['rbnd']
        o = told == C; nw = tnew == C
        drop = o & ~nw; add = (~o) & nw
        dis = (told != tnew) & ~o & ~nw          # both decided, different verdicts: impossible
        bad += int(dis.sum())
        print(f"{n:6d} {len(told):10d} {int(o.sum()):8d} {int(nw.sum()):8d} {int(drop.sum()):8d} "
              f"{int(add.sum()):7d} {(rb[drop].max() if drop.any() else 0):13.3e} "
              f"{rb.max():11.3e}" + (f"  DISAGREE {int(dis.sum())}" if dis.any() else ""))
        sys.stdout.flush()
    return bad

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--logged', action='store_true')
    ap.add_argument('--log-path', default='cusps/interval_checks.log')
    ap.add_argument('--sweep', type=int, nargs='*')
    a = ap.parse_args()
    fails = 0
    if a.logged: fails += logged(a.log_path)
    if a.sweep:  fails += sweep(a.sweep)
    sys.exit(1 if fails else 0)
