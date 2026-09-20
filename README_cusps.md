# Cusp-point pipeline

`cusps_fast.py` finds all "cusp points" (tie points p* = p*(n,i,j) with i+j>n at which the
ordered-binomial expectation E(n,p) has a local minimum) and records p*, E, F3, and the one-sided
slopes.  Definitions, method and output columns are in the script's docstring.

## Run (local machine)
    python3 -m venv .venv            # on this Mac, build from an arm64 Python (see CLAUDE.md)
    .venv/bin/pip install numpy mpmath numba
    .venv/bin/python cusps_fast.py --nmax 2000 --workers 8 --out cusps/
    .venv/bin/python cusps_fast.py --merge --out cusps/        # writes cusps/cusps_all.csv

Resumable: one file per n in cusps/; re-running skips finished n.  The ETA printed at startup is too optimistic when resuming, because the
finished n count as done work.
`cusps/interval_checks.log` lists every tie point that needed interval arithmetic and its verdict.

Merge options: `--split-mb 25` splits the output into parts (with `cusps/cusps_manifest.txt`),
`--merge-nmin/--merge-nmax` restricts the n range, `--slim` drops the S_minus/S_plus columns.

## Timings (8 workers, Apple M1 Max)
| run | time |
|---|---|
| n<=1000 | 1 min 47 s |
| n=1001..2000 (`cusps_fast.py`) | 37 min, ~40 s per n at n~2000, ~75 interval checks per n |
| further (extrapolated) | ->3000 ~3.4 h, ->4000 ~11 h, ->5000 ~26 h (probably low) |

## Current outputs
- `cusps_n2000.csv.gz`: n=3..2000, 707,417 cusps (same data as `cusps/cusps_all.csv`)
- `cusps_n1000.csv.gz`: n=3..1000, 176,816 cusps

## Recompute any row at high precision
    .venv/bin/python cusps_fast.py --recheck 34 12 30 --dps 50

## Columns of cusps_all.csv
n,i,j,pstar,E,F3,F3_sign,S_minus,S_plus,slope_left,slope_right,certified_by
(pstar/E/F3/slopes in double precision; the cusp decision itself is certified with margin 1e-6 or
interval arithmetic.  certified_by is `double` or `iv50`/`iv100`/`iv200`.)

## History
An earlier numpy-only script, `cusps_parallel.py`, produced n<=1000 but could not go further: its
near-tie test sent a tie point to interval arithmetic whenever any two adjacent sorted masses were
within relative 1e-8, and above n~1090 that happens more and more often (18% of tie points at n=2000), almost always
between deep-tail masses (~1e-300) that cannot affect the answer.  Checks went from 12 at n=1000 to 182,542 at n=2000
(~12 h per n).  `cusps_fast.py` ignores masses below `TINY = 1e-290` and ranks in O(n) without
sorting.  Rerunning n=3..1000 with it reproduced every per-n file and the interval-check log of
the old script byte for byte, so the old script was removed.
