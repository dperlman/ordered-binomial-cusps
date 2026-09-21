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
- binom_core.py is the ONLY implementation of the mathematics: constants (MARGIN/GAP/TINY), lnC,
  E_half, the numba screening kernel (tie_kernel/screen), certify()/certify_escalating(),
  evaluate(), recheck().  Every other script imports it.  Do not re-derive any of this elsewhere --
  before this existed the mass recurrence, the rank merge, the F3 formula and E_half were each
  written twice, and the two E_half versions had already drifted apart (500x different error at
  n=2000).  numpy + numba + mpmath only, so the certified numerics stay auditable.
- _one_tie() is the single place masses and ranks are computed; tie_kernel loops it and evaluate()
  calls it once.  Masses are ALWAYS normalised by their own sum (see RESEARCH_LOG section 7): E then
  matches a 60-digit computation to ~2e-13, and E - E(1/2) carries ~7.4 digits at n=2000 instead of
  ~3.6.  tie_kernel's one flag is collect_all: False returns only MIN/CHECK tie points, True returns
  every tie point.  There is deliberately NO flag for reproducing the older unnormalised output --
  flags that exist only to preserve a less accurate calculation were removed on 2026-09-20; files
  written before then differ in the last digits of E and the slope numerators, so regenerate rather
  than mixing them.
- cusps_fast.py: the certified generator (CLI, parallel driver, per-n CSV, merge, recheck).  Parallel, resumable generator of all cusp points up to
  --nmax; its docstring holds the definitions, method, column list and history.
  Run: python cusps_fast.py --nmax 2000 --workers 8 --out cusps/  then  --merge --out cusps/.
  Merge options: --split-mb, --merge-nmin/--merge-nmax, --slim.  --recheck n i j gives 50-digit values.
- numba screening kernel (recurrence for the masses, two-pointer merge for ranks, O(n) per tie
  point); masses below TINY=1e-290 are treated as zero, which keeps the interval-arithmetic
  workload at 4-86 checks per n (mean 33) for n=1001..2000.
- One CSV per n in cusps/; cusps/cusps_all.csv after merge; interval_checks.log lists every
  tie point that needed interval arithmetic and its verdict.
- dump_ties.py (binom_core + pyarrow): builds the Parquet plotting datasets for a given n --
  data/ties/n=NNNNN/ (every tie point) and data/cusps/n=NNNNN/ (cusp subset).  Run
  .venv/bin/python dump_ties.py --n 100 --data data/ [--verify].  It CERTIFIES IN PLACE using
  binom_core.certify, and records decided_by per row.  It does not look is_cusp up in
  cusps_all.csv: a lookup silently marks every tie point as a non-cusp for any n the CSV does not
  cover, which looks like a result rather than an error.  --verify cross-checks against the CSV
  where it exists, as a regression test rather than a dependency.
- Every tie point carries the full slope decomposition (A/T/V and their slopes), not just cusps;
  only the cusp-to-cusp distance columns are cusp-specific.
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

## Data publication policy (three tiers) -- keep to this
The public repo is github.com/dperlman/ordered-binomial-cusps.  Generated data is large and
reproducible, so it is NOT committed.  Anything new follows one of three tiers:

- Tier 1 -- committed, always, as CSV text (~0.5 MB total).  Small curated summaries that let
  someone reproduce the plots without downloading anything: public/per_n_summary.csv (one row per n)
  and public/cusps_decade.csv (cusps for n=100,200,500,1000,2000).  Built ONLY by
  make_public_data.py so they are reproducible and auditable, never hand-edited.  Text, not Parquet:
  line diffs keep history small when they are regenerated.
- Tier 2 -- committed, but rarely.  Mid-size derived tables (~1 MB) such as
  public/negF3_neighbors.csv.  Re-commit only when the results actually change, not after every
  rerun; each new version costs its full size in history forever.
- Tier 3 -- GitHub Releases, NEVER in git history.  The big archives (cusps_n1000.csv.gz 11 MB,
  cusps_n2000.csv.gz 45 MB) and per-n Parquet tie files.  Release assets live outside the repo, so
  they do not affect clone size, and they get stable download URLs.  Tag one release per milestone
  (v1.0-n2000 = the n<=2000 tables); a longer run means a NEW release, not a new version of a file.

Rules that keep this working:
- Never git-add anything under cusps/, data/ (except manifest.csv), analysis/, or a *.csv.gz.
  .gitignore enforces this; do not override it with `git add -f`.
- git-lfs is deliberately NOT used: on a public repo its bandwidth quota is billed to the owner for
  everyone else's clones, and it breaks clones for anyone without lfs installed.  Use a Release.
- Regenerate rather than restore.  Every ignored path has its rebuild command in .gitignore.
- If an artifact is too big for Tier 1 but someone needs it offline, it is Tier 3, not Tier 2.

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
