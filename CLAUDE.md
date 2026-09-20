# Project: cusp points of the ordered binomial distribution

READ RESEARCH_LOG.md FIRST: it holds every established result, the proof reductions, the failed
approaches and the open questions.  Append new results to it (with the n-range they were checked on).

## Background (math)
- f_p(k) = C(n,k) p^k (1-p)^(n-k), k=0..n, are the binomial masses.
- The "ordered binomial" expectation is E(n,p) = sum_k w_k f_p(k), where w_k is the rank of f_p(k)
  among all n+1 masses in increasing order (0 = smallest).  Equivalently
  E(n,p) = sum_{k<l} max(f_k, f_l) = n/2 + (1/2) sum_{k<l} |f_k - f_l|.
- E(n,p) is continuous, symmetric about p=1/2, concave on each interval between "tie points",
  with convex kinks at tie points.  All local minima are at tie points.
- A tie point p*(n,i,j), 0<i<j<n, is the p at which f_p(i)=f_p(j):
  rho = p*/(1-p*) = (C(n,i)/C(n,j))^(1/(j-i)).  We only study i+j>n (p*>1/2) by symmetry.
- A "cusp point" is a tie point that is a local minimum of E(n,.).  Test: with the left-side
  ranking (w_j = w_i - 1), S_- = sum_k w_k f(k)(k - n p*) and S_+ = S_- + (j-i) f(i);
  cusp <=> S_- < 0 < S_+.  One-sided slopes are E'_± = S_± / (p* (1-p*)).
- F3(n,i,j) = (n+i-j)(i+j-2np*) + (j-np*).  Open question: sign of F3 at cusp points.
  Known (n<=200): 7048 cusps, 70 with F3<0; all F3<0 cusps are unusually close (in p) to
  another cusp, which is always a narrow near-mode tie.
- Conjecture under study: E(n,p) >= E(n,1/2) for all p (p=1/2 is the global minimum).

## Code
- cusps_fast.py is the only script.  Parallel, resumable generator of all cusp points up to
  --nmax; its docstring holds the definitions, method, column list and history.
  Run: python cusps_fast.py --nmax 2000 --workers 8 --out cusps/  then  --merge --out cusps/.
  Merge options: --split-mb, --merge-nmin/--merge-nmax, --slim.  --recheck n i j gives 50-digit values.
- numba screening kernel (recurrence for the masses, two-pointer merge for ranks, O(n) per tie
  point); masses below TINY=1e-290 are treated as zero, which keeps the interval-arithmetic
  workload at 4-86 checks per n (mean 33) for n=1001..2000.
- One CSV per n in cusps/; cusps/cusps_all.csv after merge; interval_checks.log lists every
  tie point that needed interval arithmetic and its verdict.
- dump_ties.py (numpy+numba+pyarrow): builds the Parquet plotting datasets for a given n --
  data/ties/n=NNNNN/ (every tie point) and data/cusps/n=NNNNN/ (cusp subset).  Run
  .venv/bin/python dump_ties.py --n 100 --data data/.  It takes is_cusp from cusps_all.csv and
  never re-decides certification.
- cusps_data.py: reader for those datasets; derives width, band, w_i, f_i, T, A, V, all the slopes,
  D and the gap columns.  Use it rather than recomputing -- it encodes two numerical rules
  (normalised masses, and never subtracting S_plus-S_minus).  See RESEARCH_LOG.md section 7.
- analyze_cusps.py (numpy only): analyses straight off cusps_all.csv; writes analysis/.
- Cusp decisions are certified (double-precision screen with margin 1e-6, mpmath interval
  arithmetic for borderline cases).  Descriptive columns (E, F3, slopes) are double precision.
- Validated: n<=200 reproduces an independent 50-digit run exactly.
- Environment: .venv (Python 3.13 arm64, built from the miniconda obd env) with numpy, mpmath,
  numba.  Use .venv/bin/python.  numba caches compiled code in __pycache__; if the module is ever
  imported under a synthetic name (importlib spec_from_file_location), delete the .nbi/.nbc files
  afterwards or the next normal run fails with "No module named '<dynamic>'".

## History: why cusps_fast.py replaced the original script
- The original generator, cusps_parallel.py (numpy only), ranked the masses with a full sort and
  sent a tie point to interval arithmetic whenever ANY two adjacent sorted masses were within
  relative GAP=1e-8.  Sorting interleaves the two tails, so for large n such near-coincidences are
  unavoidable, almost always between deep-tail masses (~1e-300) that cannot move S_-/S_+.
- It produced n<=1000 in 6 min (12 checks at n=1000), but the checks switch on at n~1090 and then
  explode: 3,510 at n=1200, 23,596 at n=1400, 182,542 at n=2000 (~12 h per n, ~15 days for n<=2000).
- cusps_fast.py keeps the same criterion but ignores masses below TINY, and replaces the sort with
  an O(n) merge.  n=1001..2000 took 37 min.
- Before the old script was deleted (2026-09-18), cusps_fast.py was rerun for all n=3..1000: all 998
  per-n files and the interval-check log were byte-identical to the old script's output.

## Status
- Cusp tables complete for n<=2000 (cusps_n2000.csv.gz).  Results, counts and timings are in
  RESEARCH_LOG.md section 7.

## Conventions
- Never assume "double ties" (two pairs with the same p*) — none exist for n<=1000 (checked
  separately by the user), but flag any if found.  (For n<=2000 no two cusps share a p* at double
  precision; that is not a full check over all tie points.)
- Dependencies: install freely into the project venv (.venv) with .venv/bin/pip, within reason --
  well-known, actively maintained packages that earn their place (numpy, mpmath, numba, pyarrow,
  matplotlib and the like).  Prefer a package over hand-rolling something it already does well.
  NEVER install into the system Python or into the miniconda envs; .venv is the only target, and it
  is disposable -- it can be rebuilt from scratch at any time.  Say in your reply what you installed.
  Keep the certified generator (cusps_fast.py) lean -- numpy, mpmath, numba only -- so its numerics
  stay easy to audit; analysis, export and plotting scripts may use more.
- Keep text outputs as CSV with a header; bulk numeric data for plotting goes in Parquet (see the
  data architecture notes in RESEARCH_LOG.md).
- When adding analyses (nearest-neighbor cusp gaps, F3 sign tables, etc.), write them as
  separate small scripts that read cusps/cusps_all.csv rather than recomputing cusps.
