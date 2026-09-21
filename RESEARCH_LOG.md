# Research log: the ordered binomial distribution and its cusp points

Written 2026-09-18 from a claude.ai project conversation (Sept 9–18, 2026). Everything below was
established in that conversation; "verified" means checked numerically with the stated range,
"proved" means an argument was given and checked. Keep appending to this file.

## 1. Definitions and notation

- f_p(k) = C(n,k) p^k (1-p)^(n-k), k = 0..n.  q = 1-p.  K ~ Bin(n,p).
- Ranks: w_k = rank of f_p(k) among the n+1 masses in increasing order, 0 = smallest.
- Ordered-binomial expectation:  E(n,p) = sum_k w_k f_p(k).
- G(p) = sum_{0<=k<l<=n} |f_p(k) - f_p(l)|.   Then  G = 2E - n.
- Tie point of the pair (i,j), 0<=i<j<=n:  p*(n,i,j) = the p with f_p(i)=f_p(j):
  [Convention widened 2026-09-21 from 0<i<j<n; see the section 8 entry of that date.  The pairs
  (i,n) are real tie points with p*>1/2; the last of them, (n-1,n), sits at p*=n/(n+1), above
  which E = n p exactly.]
      rho := p*/(1-p*) = ( C(n,i)/C(n,j) )^(1/(j-i)).
  By symmetry p <-> 1-p, (i,j) <-> (n-j,n-i), only i+j>n (p*>1/2) is studied.
- F3(n,i,j) = (n+i-j)(i+j-2np*) + (j-np*).   (Question (3) of the original note, with w_i = n+i-j.)
- Cusp point = tie point at which E(n,.) has a local minimum.

## 2. Structural facts (all proved unless noted)

1. At a tie point, i+1 <= (n+1)p* <= j, with equality only for j=i+1 (then p*=(i+1)/(n+1)).
   Hence i < np* < j always.  (Answers Q1: yes.)
2. w_i = i + n - j at a tie point (all k in (i,j) have larger mass, all k outside [i,j] smaller).
3. E(n,p) = sum_{k<l} max(f_k,f_l) = n/2 + (1/2) sum_{k<l} |f_k - f_l|.  Global formula: no sorting,
   no tie convention needed.  E(n,p)=E(n,1-p).  E(n,0)=E(n,1)=n.
4. Crossing a tie point (i,j) from left to right adds exactly f_p(j)-f_p(i) to the polynomial piece.
   E(n,p) = P_0(p) + sum_t (f_p(j_t)-f_p(i_t))^+  over tie points t in (1/2,1).
5. P_0(p) (the piece on (1/2, p*_1)) = 2 E_p[min(K,n-K)] + P_p(K>n/2).  The first tie point above 1/2
   is the innermost pair: (n/2, n/2+1) for n even (p*_1=(n+2)/(2n+2)); (m, m+2) for n=2m+1
   (p*_1 = sqrt(m+2)/(sqrt m + sqrt(m+2)))  [verified n<=200].
6. Local structure at a tie point: E = B(p) + (1/2)|f_j - f_i|, with B smooth.  The V is symmetric to
   first order with arm slopes -/+ D/2, D = (j-i) f(i)/(p*q*) = the full slope jump.
   Split B = A + T, A = sum_{k not in {i,j}} w_k f_k (others), T = (r+1/2)(f_i+f_j) (pair's smooth part).
   One-sided slopes E'_± = S_±/(p*q*), S_- = sum_k w^-_k f(k)(k-np*) with w^-_j = w^-_i - 1,
   S_+ = S_- + (j-i) f(i).   Cusp  <=>  -(j-i) f(i) < S_- < 0  <=>  |B'(p*)| < D/2.
7. E is concave on every piece between tie points (verified n<=100; E'' ~ -10^3 at n=34, -10^4 at
   n=100).  So all local minima are cusps and all interior local maxima are smooth points.
   p=1/2 is a local minimum (verified n<=40) and appears to be the global minimum (the conjecture).
8. Meaning of F3:  F3 f(i)/(p*q*) is the right-hand slope of the PAIR's own contribution T+V;
   (F3-(j-i)) f(i)/(p*q*) is its left-hand slope.  Equivalently  S_+ = R + f(i) F3,
   S_- = R + f(i)(F3-(j-i)), R = contribution of the other masses.
   F3 = (2w_i+1)[x̄ - np*], x̄ = (w_i i + (w_i+1) j)/(2w_i+1) (rank-weighted position of the pair).
   Cusp condition in terms of R:  -f(i)F3 < R < f(i)(j-i-F3).
   F3<0  <=>  the pair's contribution is decreasing on both sides; a cusp then exists only because
   the other masses' curve rises steeply there (always just before a smooth local max).
9. Q2 (2np* < i+j) is FALSE in general even with p*>1/2: smallest counterexample n=9,i=2,j=8.
   It holds for adjacent ties and generally for narrow ties, fails for wide ones.
10. PROVED (externally, 2026-09-21; the user's colleagues -- ask the user for the paper): E has NO
    smooth local minima.  Every local minimum of E(n,.) is at a tie point, i.e. is a cusp.  This
    was previously listed above as a consequence of the numerically-verified concavity (fact 7);
    it now stands on its own proof and no longer depends on fact 7.
11. PROVED (here, 2026-09-21): no tie point is a local MAXIMUM -- there are no "anti-cusps".
    E = n/2 + (1/2) sum_{k<l} |f_k - f_l|.  At a tie point of the pair (i,j) the only non-smooth
    term is |f_i - f_j| = |g(p)| with g(p*) = 0 and g'(p*) != 0, and |g| has a CONVEX kink there:
    its one-sided slopes are -|g'| and +|g'|.  So E'_+ - E'_- = D = (j-i) f(i)/(p* q*) > 0 at
    every tie point (a multi-tie sums positive kinks, still > 0).  A local maximum would need
    E'_- >= 0 >= E'_+, i.e. E'_+ <= E'_-, contradicting D > 0.  Hence every local maximum of E is
    a smooth point.  With fact 10: local minima <=> tie points that are cusps; local maxima <=>
    smooth points.  No concavity needed.  (D > 0 is also verified on every stored tie point.)
    Consequence for the pipeline: nothing to add.  There is no anti-cusp to search for, no new
    column or file, and the certification criterion S_- < 0 < S_+ already captures every extremum
    that can sit at a tie point.

## 3. Numerical results on cusp points

- Pipeline: cusps_fast.py (numba kernel; O(n) per tie point via unimodal two-pointer merge).
  Certification: double-precision screen with margin 1e-6 on S_±, mpmath interval arithmetic
  (50/100/200 digits) for anything within the margin or with near-coincident other masses.
  Masses below 1e-290 are treated as zero; kink height (j-i)f(i) computed in log space.
  Measured |S_double - S_exact| ~ 6e-9 at n=2000 (grows ~ n^1.8): margin is safe to n~5000+.
- n<=200: 7048 cusps (p*>1/2).  Count per n ~ 0.355 n.  All cusp p* in (0.5025, 0.657); the max
  cusp p* is ~0.65 and roughly constant in n.
- F3 at cusps: 6978 positive, 70 negative (first at n=34: (34,12,30)).  Negative-F3 cusps sit in
  p* in (0.563,0.651), i/n in (0.32,0.54), j/n in (0.72,0.88): wide ties.  They are shallow dips
  (both slopes tiny; e.g. (34,12,30) dip depth ~3e-7, width ~6e-5 in p).
- "Double cusp" claim (adjacent tie point also a cusp): FALSE — only 12 of 70 are double cusps,
  58 isolated; same 17% rate as positive-F3 cusps.
- "Too close" claim: TRUE in the soft form.  Nearest-cusp gap for negative-F3 cusps:
  n*gap median 0.044 (max 0.161) vs 0.49 for positive-F3; intervening non-cusp ties median 4
  (max 18) vs 43; percentile within same n max 47.  No outliers.  The nearest cusp is always a
  narrow near-mode tie (e.g. (49,51), (43,47)).  Closeness is necessary but not sufficient
  (many positive-F3 cusps are equally close).  Suggested definition: gap < 1/3 of median gap at n.
- Slope jump D at cusps ~1.87 (median), nearly independent of n; one-sided slopes small.
- Closest cusp to E(1/2) is always in the first band (i+j=n+1, p* ~ 1/2 + 1/(2n)); the gap
  E(p*)-E(1/2) shrinks like ~0.04/n.  The conjecture is delicate only within O(log n / n) of 1/2.
- Tie points cluster in bands: pairs with i+j = n+m sit near p* ~ 1/2 + m/(2n) (m=1 band exactly,
  width ~ ln n/(4n)).
- No double ties (distinct pairs with identical p*) for n<=40 (exact search); user checked n<=1000.

## 4. Proof approaches for  G(p) >= G(1/2)  (equivalently E(n,p) >= E(n,1/2))

FAILED:
- Majorization / Schur-convexity: needs Top_s(p) >= Top_s(1/2) for all s (Top_s = mass of the
  s largest masses = best window of s consecutive integers).  Fails at s=1 for n=2, p=0.6.
  Fails exactly for s with parity making the optimal window at 1/2 symmetric about n/2.
- Symmetrizing f_p first: makes G smaller (n=2: 0.44 < 0.5).  The asymmetric mirror-pair terms
  |f_k - f_{n-k}| are essential.
- Pairing Top_s + Top_{s+1}: fails for n>=6 (in the first band).  Partial sums over s from
  either end: fail.
- The binomial ordering is NOT always a "distance ordering" (closer to some centre c => larger);
  fails already at n=3.

PROMISING (the main reduction):
- For any c with 2c not an integer, G(p) >= L_c(p) := sum_{k<l} sgn(2c-k-l)(f_l - f_k).
  For c in Z/2 + 1/4 the coefficient of f_k is n+1-4|k-c| (larger when the mirror 2c-k falls
  outside [0,n]).  Hence
        G(p) >= (n+1) - 4 E_p|K - c|   for every c in Z/2 + 1/4,
  with EQUALITY at p=1/2, c=n/2+1/4.
- Define D(p) = min_{c in Z/2+1/4} E_p|K-c|.  Then
        D(p) <= D(1/2) for all p   ==>   G(p) >= G(1/2) for all p.        (★)
  (★) verified on a 1200-point grid for all n<=200; it keeps 100% of the true margin for n<=20
  and ~2/3 of it at n=200 (tight spot p ~ 1/2 + 1/(2n)).
- D explicitly: D(p) = E_p|K-med| + (1/4)( P(K=med) - |P(K<med) - P(K>med)| ), med = a median.
  D(1/2) = E_{1/2}|K-n/2| + (1/4) P_{1/2}(K=n/2) [n even];  E_{1/2}|K-n/2| = (n/2) C(n-1,floor(n/2)) 2^{1-n}.
  de Moivre: E_p|K-np| = 2npq P(Bin(n-1,p) = floor(np))  (verified).
- D(p) is NOT monotone: for fixed c, phi_c(p)=E_p|K-c| has slope Cov_p(|K-c|,K)/(pq) and turns up
  before the median jumps; D is the lower envelope (saw-tooth).  So (★) <=> at every switch point
  p_m (median jumps m -> m+1), phi_{m+1/4}(p_m) <= D(1/2).  The delicate case is the first switch
  (m = n/2), slack O(1/n).  Plan: (i) epsilon-expansion of phi_{n/2+1/4} around p=1/2 for the first
  switch; (ii) asymptotics with explicit constants + finite check for later switches.

## 5. Open questions / next steps

- Prove (★), starting with the first-switch case.
- Extend cusp tables to n<=5000; check whether negative-F3 cusps remain "close" and whether max
  cusp p* stays ~0.65.  (n<=2000 done; tie-point dumps exist at n=3000,4000,5000.)
- DECIDED 2026-09-21: the tie-point convention is 0<=i<j<=n.  The kernel (cusps_fast.py,
  dump_ties.py) still implements j<n and must be updated, then the Parquet dumps rebuilt and the
  tie-point release re-issued; n=3's cusp (1,3) added.  Sequenced after the performance work so
  that work can be validated byte-for-byte first.
- Decide on the sharpened CHECK trigger (section 8): prototyped, 99% reduction with 0 disagreements
  at n<=1000, but the error constants are not yet derived and it is NOT in use.
- Conjecture (section 8, 2026-09-21): a cusp's pair mass f(i) is bounded below, ~1e-7 for n<=3000.
  How that floor moves with n is NOT established -- the low tail fits no law (see the correction
  there); only the bulk is clean, at n^-1/2.  If a bound were provable it would let the generator
  skip ~90% of tie points outright.  Rests on piecewise concavity (fact 7) and on tie-point spacing from the smooth
  maxima of E, neither of which is proved.
- Quantify cusp "depth" (dip height before the nearest smooth max) for all cusps; conjecture:
  negative-F3 cusps are the shallowest.
- Prove piecewise concavity of E (fact 7 above is numerical).
- Prove that the first tie point above 1/2 is the innermost pair (fact 5 is numerical for n<=200).

## 6. Files

- cusps_fast.py: the generator (numba); its docstring has the definitions and method.  The older numpy
  script cusps_parallel.py was removed 2026-09-18 after cusps_fast.py reproduced all of its n<=1000
  output byte for byte (see section 8).
- cusps/nNNNNN.csv: per-n results; cusps/interval_checks.log; merge with --merge [--split-mb --slim].
  cusps/cusps_all.csv = cusps_n2000.csv.gz (n<=2000); cusps_n1000.csv.gz (n<=1000).
- Columns: n,i,j,pstar,E,F3,F3_sign,S_minus,S_plus,slope_left,slope_right,certified_by
  (slim: n,i,j,pstar,E,F3,F3_sign,slope_left,slope_right).
- analyze_cusps.py (numpy only): analyses over cusps_all.csv; writes analysis/ (per_n_summary.csv,
  negF3_neighbors.csv, all_neighbors.csv, summary.txt).  Copies of per_n_summary.csv and
  negF3_neighbors.csv are at the top level.
- Original question sheet: 3QuestsClaude.pdf (page 21 of a larger note).

## 7. Data architecture (plotting datasets)

Layout, built on demand per n by dump_ties.py, read by cusps_data.py:

    data/ties/n=NNNNN/part.parquet     every tie point (p*>1/2) for that n
    data/cusps/n=NNNNN/part.parquet    the cusp subset + cusp-to-cusp distances
    data/manifest.csv                  n, n_ties, n_cusps, E_half, schema_version, seconds, built

Parquet, zstd + BYTE_STREAM_SPLIT, dictionary encoding OFF (it silently defeats BSS).  Hive
partitioning by n lets a plot scan a subset of n; column pushdown makes a two-quantity scatter read
only two columns (~0.6 ms vs 16 ms for npz).  Measured 45 B/row float64 vs 50.8 npz vs 156 CSV.
Ties and cusps share one schema and one code path, so they cannot disagree; the cusp files exist
only because "cusps across all n" would otherwise touch 16 GB of tie files to extract 0.07% of it.

Sizes: cusps for all n<=5000 ~199 MB.  Ties are per-n: 6 MB at n=1000, 24 MB at 2000, 150 MB at
5000 -- fine individually, but ALL n<=5000 would be 250 GB, so only build the n you need.

Stored (float64 unless noted): i, j (int16; n is the partition key), pstar, ln_fi, E_minus_Ehalf,
S_minus, F3, is_cusp (bool), decided_by, gap_prev, gap_next, rank_in_n (int32).
Cusp files add cusp_gap_prev/next, cusp_intervening_prev/next, nb_i, nb_j.
Derived in cusps_data.py: n, width=j-i, band=i+j-n, w_i, f_i, u=f/(p*q*), kappa, S_plus, E, T, A, V,
slope_left/right, slope_T, slope_A, slope_V_right, D, gap_nearest, n_gap_nearest, cusp_gap_nearest.

Two decisions that matter numerically, both found while testing n=100:
- Masses are normalised by their own sum before E is accumulated.  Computed as exp(lnC + k ln p +
  (n-k) ln q) they carry a shared relative error ~6e-13 (they sum to 0.999999999999363), which is
  ~1.2e-9 on E at n=2000 and leaves only ~3.6 digits on E - E(1/2) (~5.4e-6).  Normalised: ~6.3
  digits.  The summation method is irrelevant; fsum changes nothing.  E(n,1/2) uses the same
  convention, and only E_minus_Ehalf is stored (E_half is one number per n, in file metadata).
- S_plus is NOT stored.  The kink kappa = (j-i)f(i) is usually many orders below S_minus (S ~ 1.5
  with kappa ~ 1e-7, down to 1e-25), so a stored double S_plus loses it entirely: 72 of 2401 rows at
  n=100 had S_plus == S_minus bit-for-bit.  ln_fi determines kappa exactly, so the loader derives
  S_plus = S_minus + kappa and slope_right = slope_left + D, D = kappa/(p*q*).  Never compute
  slope_right as S_plus/(p*q*).  D spans 26 orders of magnitude at n=100 alone and is exact this way.

is_cusp is certified in place by dump_ties.py using binom_core.certify -- the same code the
generator runs -- and decided_by records the route ('double', 'iv50'/'iv100'/'iv200', or
'UNRESOLVED').  It is NOT looked up in cusps_all.csv: that table only covers n<=2000, and a lookup
would label every tie point of any larger n a non-cusp without raising anything.  Re-certifying
costs only the CHECK-tagged tie points (75 at n=2001: 16.8s of a 41s build).  --verify compares
against the CSV where it exists: at n=100 and n=1000 the cusp set and the certification route
agree on every row.

F3 is the right-hand slope of T+V in units of u = f/(p*q*), NOT the slope of T:
    (T+V)'_+ = F3 u,  (T+V)'_- = (F3-(j-i)) u,  T' = (F3-(j-i)/2) u,  V'_± = ±(j-i)/2 u.
Verified against the certified table to 1.6e-13 relative over 3000 rows, and the slope jump
E'_+ - E'_- = D to 6e-16.

## 8. Log entries

### 2026-09-17/18 (Claude Code): cusp tables to n=2000
- n<=2000 complete: 707,417 cusps (n<=1000: 176,816), 0 UNRESOLVED.  F3<0: 7,639 overall
  (1,920 for n<=1000; 70 for n<=200, matching section 3).  n<=200 reproduces 7048 cusps exactly.
- Interval checks: n<=1000: 1,601, all NOT.  n=1001..2000: 33,188, of which 33,169 NOT and
  19 MIN (rows with certified_by=iv50).  These 19 are the first cusps that double precision could not
  certify on its own.  5 of the 19 have F3<0 (26%, vs 1.1% overall): (1167,710,812), (1483,837,1019),
  (1642,929,1076), (1674,955,1143), (1813,1052,1199).  Plausible link: F3<0 cusps sit next to narrow
  near-mode ties (section 3), which is also what produces near-coincident masses.  Not yet analysed;
  small sample, and the selection is not independent of the mechanism.
- No two cusps share a p* at double precision for n<=2000 (not a check over all tie points).
- Why the original script stopped at ~1000: its near-tie test (any two adjacent sorted masses within
  relative 1e-8) fires on deep-tail masses ~1e-300.  Checks per n: 12 at n=1000, 533 at 1100,
  3,510 at 1200, 23,596 at 1400, 182,542 at 2000 (~12 h per n).  cusps_fast.py ignores masses below
  1e-290: 4-86 checks per n (mean 33) for n=1001..2000.
- Runtime (8 workers, Apple M1 Max): n<=1000 1 min 47 s; n=1001..2000 37 min (~40 s per n at n~2000).
  Extrapolated: to 3000 ~3.4 h, to 4000 ~11 h, to 5000 ~26 h more; the last two are probably
  underestimates since checks per n grow ~n^3.4 and each mpmath call costs O(n).
  The "n<=5000 (running)" item in section 5 has not been started here.

### 2026-09-18 (Claude Code): analyze_cusps.py on n<=2000
- Closeness of F3<0 cusps: holds for almost all, but NOT for all, beyond n=200.  n*gap median 0.046
  (F3>0: 0.499).  Percentile-at-n of the nearest-cusp gap: max 47 for n<=200 (as in section 3), but
  38 of 7,639 F3<0 cusps exceed 50 (first n=366; max 95 at (1864,1151,1280)); max gap/median 1.002.
  All 38 sit at the top edge of the cusp range, p* in (0.651, 0.653); 31 of them are the largest-p*
  cusp at their n.  Their nearest cusp is always below them, a narrow tie (width 1-4), at the
  ordinary spacing n*gap ~ 0.5.  So "F3<0 cusps are close to another cusp" fails exactly at the upper
  edge.  The nearest cusp of an F3<0 cusp has width j-i median 15 but up to 88, so "always a narrow
  near-mode tie" also needs qualifying for n<=2000.
- Double cusps: 536/7,639 F3<0 (7.0%) vs 39,302/699,776 F3>0 (5.6%).
- F3<0 fraction is flat at ~1.05-1.10% across n=3-500, 501-1000, 1001-1500, 1501-2000.
- Max cusp p*: 0.65693 overall (at n=15); for each n from 250 to 2000 it is 0.6514-0.6523.
- The 19 interval-certified cusps: 10 were escalated because a slope was within the 1e-6 margin
  (min(|S_-|,|S_+|) between 1.3e-7 and 8.3e-7, no near-tie); 9 because of a near-tie between two other
  masses (relative gap < 1e-8; both slopes >= 2.3e-3).  4 of the 5 F3<0 ones are margin cases
  (tiny slopes, e.g. (1483,837,1019): |S_-|=2.5e-5, |S_+|=1.3e-7), consistent with F3<0 cusps being
  shallow dips (section 3).  This replaces the near-tie explanation suggested in the previous entry.
  Two near-tie triggers involved masses ~5e-144 and ~2e-197, so TINY=1e-290 is still conservative.
- Closest cusp to E(1/2): in the first band (i+j = n+1) for every n<=2000, and E(p*) - E(1/2) > 0 at
  every cusp (smallest 5.4e-6 at n=2000; E in double precision).  Since the local minima of E are the
  cusps, this is numerical support for the conjecture on n<=2000.
  The gap scales as n^(-3/2), not ~0.04/n: min(E - E(1/2)) * n^(3/2) = 0.483 (n even), 0.683 (n odd) at
  n=500..2000 (fit exponent -1.499; odd/even ratio 1.412, close to sqrt 2).  Values of min*n:
  n=100: 0.048, 1000: 0.0153, 2000: 0.0108 (even).

### 2026-09-20 (Claude Code): shared core, in-place certification, public tables regenerated
- binom_core.py is now the single implementation of the kernel, certify(), E_half, lnC and
  evaluate; cusps_fast.py, dump_ties.py and analyze_cusps.py import it.  Acceptance test:
  rebuilding n=3..1000 reproduces all 998 per-n files and the interval-check log byte for byte
  (1 min 48 s), as do n=1995..2000 and their 448 check lines.
- dump_ties.py certifies in place rather than looking is_cusp up in cusps_all.csv.  n=2001 (outside
  that table) builds 999,000 ties and 706 cusps from 75 certifications in 41 s; the linear fit
  predicts 707.  Under the old lookup it would have reported zero cusps and raised nothing.
- analyze_cusps.py now uses the normalised E_half, and public/ was regenerated.  Only E_half and
  min_E_minus_Ehalf changed (max rel 7.7e-4, at large n); all counts, extremes and F3 statistics are
  unchanged.  The n^(-3/2) constants are unmoved: 0.483150 (even), 0.682404 (odd), sd ~2e-4.
- CAVEAT: this fixed only half of that quantity.  analyze_cusps reads E from cusps_all.csv, whose E
  column was computed from unnormalised masses, so min_E_minus_Ehalf is ~3.7 digits at n=2000 rather
  than ~3.6 -- not the ~6.3 the Parquet path gets by normalising both.  Cross-check of the same
  quantity: n=100 differs by 1.5e-12 between the two paths, n=1000 by 3.0e-10.  For anything needing
  more than ~4 digits on E - E(1/2), use the Parquet datasets, not cusps_all.csv.  Fixing the CSV
  would mean regenerating all of cusps/ with normalise=True (~40 min) and is not done.

### 2026-09-20 (Claude Code, later): evaluate() folded into the kernel; cusps/ regenerated
- evaluate() used to recompute the masses itself (unnormalised, ranked with lexsort) -- a second
  implementation of the same mathematics, and the reason the CSV E column was ~1.2e-9 off.  It is now
  served by _one_tie, the single place masses and ranks are computed anywhere in the project.
  E against 60-digit exact: 2.8e-14 at n=200, 0 at n=1000, 2.3e-13 at n=2000.  E - E(1/2) therefore
  carries 7.4 digits at n=2000, against 3.6 before and 3.7 after the half-fix earlier today.
- Removed the flags that existed only to reproduce the older, less accurate output: normalise on
  tie_kernel/screen and on E_half, evaluate()'s own mass computation, and the arbitrary -700 cutoff
  that zeroed the kink below ~1e-304 while ln_fi still reported it (so kappa and D disagreed for the
  deepest ties; smallest kappa at n=100 is now 5.4e-26, all positive).  collect_all remains: it
  selects which rows are emitted, not how they are computed.
- cusps/ fully regenerated (41 min 28 s, 8 workers).  Unchanged: 707,417 cusps, 7,639 with F3<0,
  0 UNRESOLVED, 34,789 interval checks, n<=200 still 7048/70, and every cusp (i,j) set, pstar, F3 and
  certification route.  Changed: E and the slope numerators, by about the size of the old error.
  per_n_summary: only min_E_minus_Ehalf, median_slope_jump and min_slope_jump moved; the n^(-3/2)
  constants are 0.483152 (even) and 0.682401 (odd), unmoved within their ~2e-4 scatter.
- HAZARD found during the rebuild: after `rm -rf cusps/`, eight stale files from the 2026-09-17 run
  reappeared as "n01993 2.csv" ... "n02000 2.csv", with their original timestamps -- a sync service
  (Dropbox is running; the folder is under ~/Documents) restoring deleted copies.  They held the OLD
  unnormalised values and matched merge's n*.csv glob, so they would have been merged in as silent
  duplicate rows.  merge() now matches ^n\d{5}\.csv$ exactly and prints anything it skips.  Check
  the merged row count against the expected total after any regeneration.

### 2026-09-20 (Claude Code, later still): faster certification and parallel tie dumps
- certify() rewritten scale-free: every decision is invariant under a common positive scale, so the
  masses are taken relative to f(i)=1 and built by the kernel's recurrence.  That removes the
  binomial coefficients entirely -- ~n-digit integers that profiling showed were 72% of the routine
  at n=4000 (computing all n+1 of them is O(n^2) in bit complexity).  Cost ~n^1.53 -> ~n^1.05.
  Validated on ALL 34,789 logged checks: 0 verdict mismatches, 0 differences in precision route.
- dump_ties.py --workers N splits the i-loop into work-balanced chunks (equal tie-point counts, not
  equal i) and distributes the certifications.  Rebuilding n=2000 in parallel gave a byte-for-byte
  identical Parquet file, metadata included.  Single-n dump times: n=2000 41.5s -> 4.9s,
  n=3000 203.8s -> 16.7s, n=4000 ~15 min (modelled) -> 43.7s, n=5000 -> 92 s.
- Tie-point dumps now exist for n = 100, 200, 500, 1000, 2000, 3000, 4000, 5000.  n=5000 has
  6,245,001 ties and 1,770 cusps (the linear fit cusps(n) = 0.3538n - 0.25 predicts 1,769).
  Check counts run BELOW the n^3.45 fit at large n: 247 at n=3000 (fit 308), 637 at n=4000 (fit 830),
  1,369 at n=5000 (fit 1,792), so the runtime projections for n>=4000 are pessimistic.
- Questions asked of the n<=2000 cusp table (707,417 cusps over 1,997 values of n):
  (a) E(cusp) > E(1/2) at EVERY cusp -- 0 failures.  Tightest margin 5.4e-6 at n=2000 against ~2e-13
      numerical error.  Since the local minima of E are the cusps, this is the conjecture holding.
  (b) The first cusp (smallest p*) has the lowest E for all but FIVE n: 6, 21, 50, 76, 125.  In each
      the second cusp is lower.  Holds for every n from 126 to 2000.
  (c) E at the cusps is NOT monotone in p*: it fails at 1,963 of 1,997 n, with 52,206 of 705,420
      adjacent cusp pairs (7.4%) stepping downward.  The decreases are spread over the whole cusp
      range (p* 0.5006..0.6520, median 0.6057), not confined to an edge.

### 2026-09-20 (Claude Code): feasibility of a sharpened CHECK trigger -- NOT IMPLEMENTED
Prototype only (numpy, outside the kernel); nothing in the pipeline uses this yet.

The CHECK trigger has two parts.  The margin test (|S_+-| within MARGIN=1e-6) is honest.  The
near-tie test (two adjacent masses within relative GAP=1e-8) is a PROXY for "the ranking may be
wrong", and it never asks whether a wrong ranking would matter.  Swapping adjacent ranks of masses
a,b changes S_- by exactly f_a(a-np*) - f_b(b-np*) ~ f*(a-b): the damage scales with the SIZE of the
masses, not with how close they are.  A near-tie between two masses at 1e-144 (which really does
trigger certification today) perturbs S_- by ~1e-141 against an S_- of order 1e-2.

Sharpened criterion: flag only when a computed bound on the total possible perturbation exceeds the
distance of S_-/S_+ from zero.  Four terms, all per tie point: (1) summation rounding
(n+1)*eps*sum|w_k f_k a_k|; (2) mass error delta_f*sum|w_k f_k a_k|; (3) error in p*, into both f_k
and a_k; (4) re-ranking -- masses whose neighbouring order double precision cannot resolve are
grouped into clusters, and within a cluster of size c any rank moves by at most c-1, giving
(c-1)*sum_{k in C} f_k |a_k|.  This replaces BOTH MARGIN and GAP with computed quantities.

Measured on the 1,601 logged checks with n<=1000: 16 still flagged, 1,585 decided in double
precision, and ALL 1,585 agree with the certified verdict.  A 99.0% reduction with 0 disagreements.
At n=1500 and n=2000, every single CHECK came from the near-tie proxy -- none from the margin.

Three caveats before this can be trusted:
- delta_f = 1e-11 and delta_p = 1e-12 are ASSERTED, conservative against a measured ~1e-15, not
  derived.  Making them rigorous is the real remaining work.
- Only tie points the old rule flagged were tested.  The other direction is untested: the new bound
  CAN flag things the old rule decided (a near-mode cluster with f ~ 1e-2 and |a| ~ n gives a bound
  ~10), so the check count could rise somewhere.  Needs measuring over all tie points, which the
  numpy prototype is too slow to do (0.4 s per tie point).
- Zero disagreements is necessary, not sufficient; correctness rests on the error analysis above.

If it holds, the certification term largely disappears and a full run to n<=5000 becomes screening-
dominated, roughly 15 h instead of 3.2 days.  Validation is cheap and direct: for every tie point the
new trigger drops that the old rule flagged, certify anyway and confirm the verdict is unchanged --
34,789 such cases exist for n<=2000.

### 2026-09-21 (Claude Code): cusp tables extended to n<=3000
- Run finished in 2 h 25 m on 8 workers (predicted ~3 h; the faster certify() accounts for the gap).
- n<=3000 complete: 1,591,532 cusps, 17,081 with F3<0 (1.07%), 0 UNRESOLVED.  The new range
  n=2001..3000 contributed 884,115 cusps and 9,442 with F3<0, so the F3<0 rate is flat.
- The n<=2000 portion is unchanged, byte-identical to the archived version.
- Interval checks: 195,243 total, of which 160,454 are new (n=2001..3000) -- about 160 per n against
  77 at n=2000, following the n^3.45 growth.  Of those 160,454, only 25 were cusps: 99.98% of the
  mpmath work now goes into confirming NON-cusps, which is exactly what the sharpened trigger in the
  previous entry would remove.
- No stale duplicate files appeared this time, confirming that moving the project out of ~/Documents
  fixed the iCloud restore hazard that nearly corrupted the previous merge.
- analyze_cusps rerun: 44 interval-certified cusps overall (9 with F3<0); the closest cusp to E(1/2)
  is still always in the first band (i+j=n+1) for every n<=3000, and min (E-E(1/2))*n = 0.0088.

### 2026-09-21 (Claude Code): the p=1/2 axis is now a row in the tie-point dumps
- p=1/2 is NOT an ordinary tie point: all mirror pairs (i,n-i) tie there simultaneously (n/2 of
  them), which is exactly why the i+j>n filter excludes it -- it IS the symmetry axis.  The
  single-pair bookkeeping S_+ = S_- + (j-i)f(i) does not describe it.
- What does apply is E(p) = E(1-p), giving E'(1/2-) = -E'(1/2+) exactly, hence S_- = -S_+ and
  kink position S_-/kappa = -1/2 EXACTLY: the zero sits dead centre in the slope jump, so p=1/2 is
  the most robust cusp there is.  Verified against finite differences.
- p=1/2 is a cusp (S_- < 0 < S_+) for EVERY n from 3 to 3000, extending section 2 fact 7, which had
  it verified only to n<=40.  Smallest S_+ is 0.375 at n=3, so double precision settles it with an
  enormous margin; decided_by is recorded as 'symmetry', not a numerical route.
- CAUTION for anyone recomputing this: f(k) and f(n-k) are equal in exact arithmetic but differ in
  the last ulp through lgamma, so a stable sort ranks them by numerical noise and returns the WRONG
  SIGN for S_+.  It did, for 26 values of n, until the masses were symmetrised (f + f[::-1])/2 and
  the ranking built with lexsort((k, f)) -- ties to the smaller index, which is the ordering just
  to the right of 1/2.
- Schema version 3: each partition now carries the axis as row 0, with sentinel (i,j) = (0,n) so
  band = i+j-n = 0 marks it (every real tie point has band >= 1), plus an n_tied_pairs column.
  E_minus_Ehalf is exactly 0 there and F3 is NaN; cusps_data blanks T/A/V/w_i/f_i on that row,
  since the pair decomposition does not apply, while kappa, D and the slopes remain correct.
  All 11 partitions rebuilt (18 min): n = 100..8000, each gaining exactly one row and one cusp.
- cusps_data: the old derived name "u" meant f/(p*q*) while the plotting code used u for S_-/kappa.
  Renamed to slope_unit, and kink_pos = S_-/kappa added as a first-class derived column.

### 2026-09-21 (Claude Code): cusps need a non-negligible pair mass -- an unproved lower bound
- Mechanism (the sawtooth): on each piece between tie points E' decreases (concavity, fact 7 --
  itself only VERIFIED n<=100, not proved); at every tie point E' jumps UP by D = kappa/(p*q*),
  kappa = (j-i)f(i).  A cusp is a jump that carries E' from negative to positive, so it needs
  |S_-| < kappa: a tiny pair mass gives a tiny kick, the tie point sits on the curve and barely
  bends it.  That is why the curved sections carry thousands of tie points and few cusps.
- Observed floor over ALL cusps n<=3000: min f(i) = 6.56e-8 (n=1075, pair (572,734)), with near-
  repeats ~7e-8 at n=298 and n=2775.
  CORRECTION: an earlier statement in this session that "cusps need f(i) >~ 1e-5" came from three
  sampled n with high floors and is wrong by ~2.5 orders; 1e-7 is the right scale.
  SECOND CORRECTION (2026-09-21, plotting): the claim here that the per-n floor follows ~4.6*n^-1.67
  was overstated.  That fit has log-log correlation only -0.55, and other order statistics of the
  same low tail give incompatible exponents -- n^-0.97 for the 1st percentile (corr -0.27) and
  n^-3.32 for the 5th smallest (corr -0.78).  The low tail follows no law; it is a rare-event
  statistic.  What IS clean: the median and the maximum both scale as n^-0.50 with log-log
  correlation -1.00 (median ~ 0.593 n^-0.50), which is just the 1/sqrt(n) of the masses themselves.
  The running minimum steps down only 10 times over n=300..3000 and has been flat at 6.56e-8 since
  n=1075.  Spread within a single n is large (429x at n=3000: 3.6e-5 to 1.5e-2).
- Prize if it could be made rigorous: at n=8000, 93.5% of tie points have f(i) < 1e-5, 89.3% have
  f(i) < 1e-10, 84.1% have f(i) < 1e-20, and none of those is a cusp.  Skipping them would remove
  most of the O(n^2) tie-point evaluations -- a bigger lever than the TINY window (which only
  shortens each evaluation) or the sharpened trigger (which only cuts certification).
- Empirical margin of a skip rule, ratio = |S_-|/kappa (cusp needs < 1), over tie points below the
  threshold:  threshold 1e-10 -> min ratio 1e6 (n=1000), 1e5 (n=3000), 1e4 (n=8000);
  threshold 1e-20 -> min ratio 1e14 / 1e15 / 1e14.  The margin ERODES with n at a fixed threshold,
  so any bound must be n-dependent -- a fixed number is not safe indefinitely.
- Why no rigorous bound falls out: a tiny-kappa cusp needs the tie point within ~kappa/|E''| of a
  smooth local maximum of E.  Nothing forbids that; ruling it out needs (a) a lower bound on |E''|
  on each piece -- unproved, fact 7 is numerical -- and (b) control of how close tie points p*(n,i,j)
  can sit to the zeros of E', a spacing/Diophantine question with no theorem behind it.  So this is
  a RESEARCH item (an open conjecture worth stating), not an engineering one.  Not implemented.

### 2026-09-21 (Claude Code): slope plot at n=1000 (plotting/slope_plot.py)
- E'_- and E'_+ at all 249,001 tie points, two dots per p*, no lines, coloured by log10 f(i).
  Both x = p* and x = index versions; 10000x6000.
- What it shows: E' at tie points fills a band that fans out from the antisymmetric pair
  E' = -/+25.23 at the p=1/2 axis, spans roughly [-25, +25] near p=1/2, and lifts entirely above
  zero at p ~ 0.65 -- which is exactly where the cusps stop (max cusp p* = 0.6519 at n=1000).  Past
  that no kick can reach zero.  Toward p -> 1, E' grows to ~640 because p*q* -> 0.
- For 83% of tie points the two dots are within one pixel (kink below resolution); the visible gaps
  are the kinks, and the cusps are the 352 pairs whose gap straddles zero.
- Pair mass vs width at n=1000: corr(ln f, width) = -0.947, corr(ln f, -width^2/(8npq)) = +0.998,
  i.e. ln f ~ c - (j-i)^2/(8 n p* q*): "nearby swap = larger mass" holds with width measured in
  standard deviations sqrt(npq), which shrink toward p=1.  corr(ln f, ln kappa) = 1.000.

### 2026-09-21 (Claude Code): cusp mass floor vs n (plotting/mass_floor.py)
- All 1,591,532 cusps of n<=3000 plotted as f(i) against n, log-log, with per-n median, 1st
  percentile, minimum and the running lower envelope; plus minima and 1st percentiles at
  n = 4000..8000 from the Parquet dumps.  f(i) is recovered from the CSV as (S_+ - S_-)/(j-i);
  cross-checked against the stored ln_fi at n=1000/2000/3000, max relative difference 7e-16.
- The bulk is clean: median and max both ~ n^-1/2, log-log corr -1.00.  The low tail is not a
  trend -- see the correction in the previous entry.
- Visible structure worth following up: the cusps do not fill the band between the median and the
  floor.  They form distinct downward streaks, each spanning ~2 orders of magnitude in f(i) and
  starting at a sharply defined upper edge.  The streaks repeat across n rather than scattering,
  which suggests they are the (i,j) bands (i+j = n+m) seen in the mass coordinate.  Not analysed.
- The n>=4000 points sit consistently ABOVE the n<=3000 envelope (minima 2.4e-6 .. 2.8e-5 against
  6.6e-8 at n=1075).  That is expected from sampling only 5 values of n rather than 2,997: the
  envelope is built from rare events, so it needs many n to dip.  It is NOT evidence that the floor
  rises with n.

### 2026-09-21 (Claude Code): E at the first three cusps (plotting/first_cusps.py)
- Raw E is useless for this: E ~ n/2, and at n=1000 E(1/2) and the first three cusps all read
  975.2624.  The content is in E - E(1/2), which for the k-th cusp is a clean power law.
- Each of the first three cusps has (E - E(1/2)) * n^(3/2) CONSTANT over n=500..3000:
      1st cusp: 0.48324 (n even, sd 1.8e-4) / 0.68254 (n odd, sd 2.7e-4)   ratio 1.4125 ~ sqrt 2
      2nd cusp: 1.96298 (even, sd 8.2e-4)   / 1.96311 (odd, sd 7.5e-4)     -- PARITY-INDEPENDENT
      3rd cusp: 4.03423 (even, sd 1.2e-1)   / 4.23822 (odd, sd 6.4e-2)
  So the n^(-3/2) law found earlier for the closest cusp extends to the 2nd and 3rd, and the
  even/odd split is a property of the FIRST cusp alone: the 2nd has the same constant for both
  parities, to 4 decimal places.  Ratios to the 1st: 4.062 / 8.348 (even), 2.876 / 6.210 (odd) --
  the parity dependence of the ratios is entirely inherited from the denominator.
- The k-th cusp is the k-th band (i+j = n+k) essentially always: 100% for k=1 and k=2 over
  n=500..3000, and 99.8% for k=3.  The 5 exceptions (n = 792, 841, 890, 994, 1600) are n where the
  m=2 band contributes TWO cusps, so the third cusp by p* is still in band 2.  Those are the
  visible dips in panel (c) and the only source of scatter in the 3rd-cusp constant (sd 1.2e-1
  against 8e-4 for the 2nd).
- Not explained: why the 2nd cusp's constant is parity-free while the 1st splits by sqrt 2, and
  what 1.96298 and 4.034 are in closed form.  The n^(-3/2) scaling itself also has no proof yet.

### 2026-09-21 (Claude Code): the 0<i<j<n convention EXCLUDES real tie points (i,n)
Raised by the user asking why E/n was not the natural quantity to plot.  Their point: above the
last tie point the ranking is the natural order w_k = k, so E = sum_k k f_k = E[K] = n p, hence
E/n = p exactly and every n collapses onto the diagonal.  VERIFIED to machine precision (E - np =
0 or ~1e-14 at n=10,50,200).  But finding "the last tie point" exposed a gap.

- At n=10 the sorted order of the masses changes at 25 values of p in (0.5,1).  Our tie list
  (0<i<j<n) contains 16 of them.  The 9 missing are exactly the pairs (i, n), i = 1..n-1.
- Those are genuine order changes and genuine kinks in E.  The convention excludes i=0 and j=n; by
  the symmetry (i,j) <-> (n-j,n-i) the pairs (0,j) map to (i,n), and that whole family maps to
  itself -- but (0,j) has p* < 1/2 while (i,n) has p* > 1/2, so restricting to p*>1/2 does NOT
  dispose of them.  Excluding both loses n-1 tie points per n from our half of the domain.
- The true last tie point is (n-1,n) at p* = n/(n+1), not (n-2,n-1) at (n-1)/(n+1) as stored.
  Above n/(n+1) every mass is in natural order and E = n p exactly.
- CUSP IMPACT: none for n >= 4, now EXHAUSTIVELY certified.  All 4,498,497 excluded pairs (i,n)
  for n=4..3000 were run through the pipeline's own criterion -- double screen with MARGIN,
  escalation on margin OR near-tie, mpmath interval verdict -- in 8 min 19 s on 8 workers:
  198 needed interval arithmetic, 0 unresolved, 0 cusps.  (An earlier sampled check over
  88,493 pairs had found the same.)  The only cusp in the whole family is n=3, pair (1,3),
  p*=0.6340, confirmed a local minimum by direct evaluation of E; our table has no n=3 rows at
  all.  So cusps_all.csv is COMPLETE for n>=4 under the widened convention without regeneration.
- WHY they are not cusps, structurally: at the tie (i,n) the common mass IS f(n) = p*^n.  For
  p* > 0.66 there are no cusps at all (E' > 0 there); for p* < 0.66, max f over such ties is
  2.8e-6 at n=30, 4.8e-10 at n=50, 1.7e-181 at n=1000 -- below the observed cusp mass floor ~1e-7
  once n > ~25.  n=3 is the exception precisely because p^3 = 0.25 is not small.
- DATA IMPACT: the Parquet tie-point dumps are missing n-1 rows per n (999 at n=1000, 2999 at
  n=3000), of which 216 and 647 respectively lie inside the cusp range p*<0.66.  Any plot of "all
  tie points" is incomplete on the right-hand side, and the E/n plot needs them to reach n/(n+1).
- NOT CHANGED.  0<i<j<n is the convention in CLAUDE.md and in the original note; whether to widen
  it to 0<=i<j<=n is the user's call, not a bug to fix silently.
