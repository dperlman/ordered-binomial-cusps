# Project: cusp points of the ordered binomial distribution

READ RESEARCH_LOG.md FIRST: it holds every established result, the proof reductions, the failed
approaches and the open questions.  Append new results to it (with the n-range they were checked on).
FACTS.md is the one-page list of what is PROVEN vs SCREENED, with ranges, programs and dates.  Keep it
in step with the log: when a result is locked in or corrected, update both.  No proofs in FACTS.md.

## Background (math)
- f_p(k) = C(n,k) p^k (1-p)^(n-k), k=0..n, are the binomial masses.
- The "ordered binomial" expectation is E(n,p) = sum_k w_k f_p(k), where w_k is the rank of f_p(k)
  among all n+1 masses in increasing order (0 = smallest).  Equivalently
  E(n,p) = sum_{k<l} max(f_k, f_l) = n/2 + (1/2) sum_{k<l} |f_k - f_l|.
- E(n,p) is continuous, symmetric about p=1/2, concave on each interval between "tie points"
  (numerical, n<=100), with convex kinks at tie points.  PROVED: all local minima are at tie points
  (external paper, via the user) and no tie point is a local maximum (the kink is convex, so
  E'_+ - E'_- = D > 0; RESEARCH_LOG.md section 2, facts 10-11).  So there are no "anti-cusps" to
  look for.  NOT proved: that there are no cusps above p ~ 0.66 (observed max cusp p* = 0.657).
- A tie point p*(n,i,j), 0<=i<j<=n, is the p at which f_p(i)=f_p(j):
  rho = p*/(1-p*) = (C(n,i)/C(n,j))^(1/(j-i)).  We only study i+j>n (p*>1/2) by symmetry.
  CONVENTION CHANGED 2026-09-21: it used to read 0<i<j<n, which silently dropped the pairs (i,n),
  i=1..n-1 -- real order changes with p*>1/2 (their mirrors (0,j) sit below 1/2, so the symmetry
  restriction does not remove them).  The true last tie point is (n-1,n) at p*=n/(n+1); above it
  every mass is in natural order and E = n p exactly, so E/n = p.  Those pairs hold no cusps for
  n>=4 (the common mass there is p*^n, far below the cusp mass floor), so cusps_all.csv is sound
  for n>=4; n=3 has one cusp, pair (1,3), that the old range never produced.
  DONE 2026-09-21: the kernel implements j<=n, cusps_all.csv carries n=3's cusp, and the Parquet
  dumps were rebuilt.
- A "cusp point" is a tie point that is a local minimum of E(n,.).  Test: with the left-side
  ranking (w_j = w_i - 1), S_- = sum_k w_k f(k)(k - n p*) and S_+ = S_- + (j-i) f(i);
  cusp <=> S_- < 0 < S_+.  One-sided slopes are E'_± = S_± / (p* (1-p*)).
- F3(n,i,j) = (n+i-j)(i+j-2np*) + (j-np*).  Open question: sign of F3 at cusp points.
  Known (n<=200): 7048 cusps, 70 with F3<0; all F3<0 cusps are unusually close (in p) to
  another cusp, which is always a narrow near-mode tie.
- Conjecture under study: E(n,p) >= E(n,1/2) for all p (p=1/2 is the global minimum).

## Code
- certify() is scale-free: masses relative to f(i)=1, built by the same recurrence as the kernel,
  so no binomial coefficients (they were ~n-digit integers costing 72% of the routine at n=4000).
  Cost went from ~n^1.53 to ~n^1.05; verified identical verdict AND identical precision route on all
  34,789 logged checks.  S_+ = S_- + (j-i) exactly in these units.
- binom_core.py is the ONLY implementation of the mathematics: constants (MARGIN/GAP/TINY), lnC,
  E_half, the numba screening kernel (tie_kernel/screen), certify()/certify_escalating(),
  evaluate(), recheck().  Every other script imports it.  Do not re-derive any of this elsewhere --
  before this existed the mass recurrence, the rank merge, the F3 formula and E_half were each
  written twice, and the two E_half versions had already drifted apart (500x different error at
  n=2000).  numpy + numba + mpmath only, so the certified numerics stay auditable.
- _one_tie() is the single place masses and ranks are computed; tie_kernel loops it and evaluate()
  calls it once.  Masses are ALWAYS normalised by their own sum (see RESEARCH_LOG section 7, data architecture): E then
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
- SHARPENED CHECK TRIGGER (2026-09-21) -- THE DEFAULT.  core.screen(..., sharp=True) replaces the
  GAP=1e-8 near-tie proxy with a derived re-ranking cluster bound added to MARGIN.  It leaves every
  verdict unchanged -- validated against all 195,243 recorded interval verdicts for n<=3000
  (0 disagreements) and over ~14.6M tie points in the other direction (0 newly flagged) -- and takes
  the checks at n=5000 from 1411 to 1.  GAP is now used by nothing but the comparison path:
  cusps_fast.py --legacy-trigger restores the old proxy for reproducing pre-2026-09-21 runs, and
  gives byte-identical cusp tables, just far more mpmath.  Re-run validate_trigger.py after ANY
  change to the kernel's arithmetic.  See RESEARCH_LOG.md, the 2026-09-21 entry, for the bound, the
  one genuine re-ranking flag it found, and its caveats.
- validate_trigger.py: the regression test for that trigger.  --logged replays every escalation in
  cusps/interval_checks.log against the recorded mpmath verdict; --sweep N... compares both triggers
  over every tie point of each N.
- One CSV per n in cusps/; cusps/cusps_all.csv after merge; interval_checks.log lists every
  tie point that needed interval arithmetic and its verdict.
- dump_ties.py (binom_core + pyarrow): builds the Parquet plotting datasets for a given n --
  data/ties/n=NNNNN/ (every tie point) and data/cusps/n=NNNNN/ (cusp subset).  Run
  .venv/bin/python dump_ties.py --n 100 --data data/ [--verify] [--workers 8].  It CERTIFIES IN
  PLACE using binom_core.certify, and records decided_by per row.  --workers splits the i-loop into
  work-balanced chunks and distributes the certifications; output is identical however it is cut
  (verified byte-for-byte at n=2000).  n=5000 takes 92 s this way.  It does not look is_cusp up in
  cusps_all.csv: a lookup silently marks every tie point as a non-cusp for any n the CSV does not
  cover, which looks like a result rather than an error.  --verify cross-checks against the CSV
  where it exists, as a regression test rather than a dependency.
- Every tie point carries the full slope decomposition (A/T/V and their slopes), not just cusps;
  only the cusp-to-cusp distance columns are cusp-specific.
- cusps_data.py: reader for those datasets; derives width, band, w_i, f_i, T, A, V, all the slopes,
  D and the gap columns.  Use it rather than recomputing -- it encodes two numerical rules
  (normalised masses, and never subtracting S_plus-S_minus).  See RESEARCH_LOG.md section 7 (data architecture).
- analyze_cusps.py (numpy only): analyses straight off cusps_all.csv; writes analysis/.
- check_collisions.py / screen_collisions.py / verify_fact_c.py: the tie-point collision tools; see
  Conventions below.
- ALWAYS pass --workers to long compute.  Jobs the agent launches inherit a reduced QoS and a
  single-threaded one lands on the efficiency cores while the performance cores idle; the same
  n<=1500 scan is 15.7 s serial and 3.0 s on 8 workers (711% CPU), because the pool spills onto the
  performance cores.
- plotting/: plots to plots/.  plotting/_style.py holds the shared defaults -- DEFAULT RESOLUTION IS
  6000x3000 (changed 2026-09-22 from 10000x6000, which was larger than anything needed).  It also
  holds marker_size(), which sizes markers to ~3x the point spacing: these plots routinely put ~5000
  points on the x axis, where a 1 px marker at 1.1 px spacing produces moire.  Markers slightly WIDER
  than the spacing, alpha < 1, antialiasing on, and NEVER a connecting line through dense points.
  Default x axis is LINEAR so all ~5000 values of n get their own column of pixels; --logx for the
  power laws.
- ONE PIXEL COLUMN PER n (2026-09-23) -- the default for EVERY plot with n on a linear x axis.  Use
  plotting/_style.py NGrid: every n gets exactly k whole pixel columns (k = largest integer with
  N*k <= 6000, so 1 px per n once N > 3000), so every per-n mark is exactly k px wide.  The FIGURE
  WIDTH FOLLOWS THE n-RANGE (margins + N*k); height stays at the default.  The data are painted
  into an RGBA array and placed with figimage -- no matplotlib resampling -- via g.points() (per-n
  marks, k px wide by h px tall) and g.density() (clouds of many points per n, shaded by count per
  pixel); curves, annotations and the legend go on g.ax.  Never bbox_inches="tight"; never change
  xlim.  marker_size() is for --logx only, where uniform columns are impossible.  n ranges start at
  the first n with data -- n=3 for cusps (n=2 has a tie point but no cusp).  Look at the PNG at
  100%: a fit-to-window viewer resamples it and brings moire back.  Older scripts are NOT yet
  converted; convert one when the user next asks for that plot.
- plotting/mass_floor_linear.py: the first NGrid plot -- cusp pair mass vs n, n=3..5000.
- plotting/lowest_cusp.py: the minimum-E cusp of each n -- its E-E(1/2), its p*, and its width.
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
- Cusp tables complete for n<=5000 (cusps_n5000.csv.gz, 292 MB): 4,421,154 cusps, 47,299 with F3<0
  (1.07%), 0 UNRESOLVED.  Complete under the widened convention, n=3's cusp included.  The n<=3000
  portion is byte-identical to the archived cusps_n3000.csv.gz.  n=3001..5000 was the first range
  generated with the sharpened CHECK trigger as the default: 399 interval checks over those 2000
  values of n, and only 94 rows in the whole 4.4M table were certified by mpmath rather than double.
  Run: 9 h 08 m on 8 workers (2026-09-21/22).  cusps/ is now 1.3 GB (672 MB of per-n CSVs plus the
  672 MB merged file).
- Tie-point Parquet dumps rebuilt under the widened range (n = 100..1000 by hundreds, 2000..8000
  by thousands).  The kernel screens only the TINY window per tie point, which is byte-identical to
  the full pass and 1.1-1.6x faster over n=1000..3000.
- Tie-point Parquet dumps exist for n = 100, 200, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000
  (schema v3: row 0 of each is the p=1/2 symmetry axis; see RESEARCH_LOG.md section 8).
- Results, counts and timings are in RESEARCH_LOG.md section 8 (log entries).

## Data publication policy -- keep to this
The public repo is github.com/dperlman/ordered-binomial-cusps.  It offers exactly two things: the
CODE AND THE ACCUMULATED KNOWLEDGE to generate these tables at any n, and a FIXED-SIZE SAMPLE of
what the output looks like.  It deliberately does NOT offer the full tables.  They are gigabytes,
they regenerate from the code in hours, and when they WERE published as GitHub release assets
(5 releases, 762 MB, 2026-09-20..22) every single asset recorded ZERO downloads.  Those releases
were deleted on 2026-09-22 and the release mechanism is not used.  The local frontier (cusps/,
data/, analysis/, the *.csv.gz archives, currently ~4 GB) simply stays local.

Committed data is divided by SHAPE, not by size, because shape is what determines growth:

- RESULT -- one row per n.  Aggregates, counts, extremes; no raw rows.  Grows LINEARLY in n
  (~164 bytes per n), so it MAY track the local frontier freely: it is the finding itself.
  Currently: public/per_n_summary.csv (n=3..5000).
- SAMPLE -- one row per cusp.  Raw rows, so it shows what the output actually looks like.  Grows
  as n^2, or per n added to the list, so it is PINNED and does NOT move when the frontier moves.
  Currently: public/cusps_decade.csv (n=100,200,500,1000,2000,3000,4000,5000) and
  public/cusps_F3_negative.csv (F3<0 cusps with nearest-cusp metrics, n<=3000).

THE PINS ARE FINAL.  DECADE and F3_NMAX in make_public_data.py were set on 2026-09-22 and are not
to be raised as a routine consequence of a longer run -- that is exactly the drift the policy
exists to stop.  Raising one is a deliberate decision by the user, re-measured against the bands.

SIZE BANDS, enforced in make_public_data.py rather than merely written here (prose is advice an
agent can rationalise past; sys.exit is not):
      <= 1 MB   fine
      1-2 MB    fine, mention it in the commit message
      2-5 MB    HUMAN DISCRETION -- refused unless a PERSON passes --allow-large
      >  5 MB   HARD STOP -- does not go in git, at any size, for any reason

Rules that keep this working:
- Never git-add anything under cusps/, data/, analysis/, or a *.csv.gz.  .gitignore enforces this;
  do not override it with `git add -f`.
- public/ is built ONLY by make_public_data.py, never hand-edited, so it is reproducible.  It
  changes only for a bug fix or a column-schema change -- NOT because n grew.
- ALL committed text is LF.  Python's csv module defaults to the 'excel' dialect (CRLF) on every
  platform, so every writer passes lineterminator="\n" explicitly.  Never let CRLF back in.
- git-lfs is deliberately NOT used: on a public repo its bandwidth quota is billed to the owner for
  everyone else's clones, and it breaks clones for anyone without lfs installed.
- Regenerate rather than restore.  Every ignored path has its rebuild command in .gitignore.
- Cost check: the whole repo with full history is ~2.3 MB to clone.  Note that git stores blobs
  zlib- AND delta-compressed, so an APPENDED CSV costs far less than its raw size -- 12 versions of
  the three public tables, 10.65 MB of raw content, pack into under 2 MB.  Do not reason about
  history cost from raw file sizes; measure it with a probe clone and `git gc`.

## Conventions
- TIE-POINT COLLISIONS ("double ties": two different (i,j) at the same n with the same p*): NONE for
  3 <= n <= 100,000.  See FACTS.md S1-S3 and RESEARCH_LOG 2026-09-24 for what each range rests on.
  n <= 10,000 is confirmed by three independent methods; 10,000 < n <= 100,000 by screen_collisions.py,
  resting on Facts A, B, C (all proved), the implementation, and lgamma behaving as documented.
  Double precision CANNOT decide this question by itself: distinct tie points are bit-identical in
  double at n=7500 and n=8333.  Every float step is either a screen followed by an exact check, or a
  binary search with a provable margin.  Flag any collision found above n=100,000.
- Collision tools.  check_collisions.py --exhaustive: float screen at 1e-9 + exact p-adic check of
  EVERY pair within tol (was neighbours-only until 2026-09-24 -- unsound, fixed).  ~N^3; 26 min to
  10,000.  screen_collisions.py: the Fact C screen; --save writes the catalogue of reducing tie
  points; --no-fact-c finds reducing points by brute force instead (float-free, independent of Fact
  C, ~34 min to 10,000); --verify compares the window against the unrestricted scan (O(n^2) memory,
  small n only).  verify_fact_c.py: independent O(n)-memory brute force for spot checks at large n,
  compared against the screen and a saved catalogue.  9.4 h to 100,000 on 8 workers.
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
