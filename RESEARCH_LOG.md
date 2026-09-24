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
- **COME BACK TO THIS** (flagged by the user 2026-09-22, deferred until the data-tier housekeeping
  is finished): interpret the n<=5000 analysis.  The run is DONE and the numbers are in the
  2026-09-22 entry -- max cusp p* flat at 0.6521-0.6523 from n=1250 to 5000 with no drift toward
  0.66; F3<0 rate pinned at 1.07% across every scale; the F3<0 "always close to another cusp"
  signature unchanged over 2.8M new cusps (median n*gap 0.045 vs 0.499, median 101 intervening ties
  vs 1109); closest cusp to E(1/2) still always in band i+j=n+1 and min (E-E(1/2))*n down to 0.0068.
  These are interesting and have NOT been discussed or followed up.
- DONE 2026-09-22: cusp tables extended to n<=5000; negative-F3 cusps do remain "close" and max
  cusp p* does stay ~0.652.  (See the entry for that date.)
- DECIDED 2026-09-21: the tie-point convention is 0<=i<j<=n.  The kernel (cusps_fast.py,
  dump_ties.py) still implements j<n and must be updated, then the Parquet dumps rebuilt and the
  tie-point release re-issued; n=3's cusp (1,3) added.  Sequenced after the performance work so
  that work can be validated byte-for-byte first.
- DONE 2026-09-21: the sharpened CHECK trigger is derived, validated (0 disagreements on all
  195,243 recorded verdicts; 0 newly flagged over ~14.6M tie points) and is now the DEFAULT.
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

### 2026-09-21 (Claude Code): slope-jump plots; the OBD repo's slope data are suspect
- plotting/slope_jump.py: D = E'_+ - E'_- at every tie point vs tie-point index, coloured by pair
  mass, log and linear variants.  D > 0 at every tie point (fact 11); verified independently by
  building the right-side ranking (w_j = w_i + 1) directly and subtracting -- 0 negatives at n=100
  and n=1000 -- since in this pipeline S_+ is derived as S_- + kappa and cannot show one.
- At n=100 on a linear axis: a sharp upper envelope at D ~ 2 that the cusps ride (D ~ 1.9 early,
  drifting to ~0.7 at the end of the cusp range), rising to D ~ 18 as p -> 1 because of the
  1/(p*q*) factor; the p=1/2 axis sits at D = 2E'_+ = 16.  Everything else is compressed at 0.
- FLAG for the OTHER repository, ~/git/OBD (OBDsaveSourceData.py and whatever consumes it): its
  plots of "right slope minus left slope at tie points" show NEGATIVE values.  That is impossible
  (fact 11), so those slope datasets and plots are wrong and should be regenerated from this
  pipeline or fixed.  Not yet diagnosed by reading the code.  Leading hypothesis: slopes estimated
  by finite differences of E, where any stencil not bracketing p* exactly measures the concave
  curvature of the smooth piece (~E'' dp < 0) instead of the kink, which dominates whenever the
  true D is small -- i.e. the negatives should sit on low-mass tie points and never on cusps.
  If instead the negatives appear on the envelope or at cusps, the old code computes something
  else and needs reading.

### 2026-09-21 (Claude Code): the cusps with a tiny slope jump ARE the F3<0 cusps
- In the linear slope-jump plots (plotting/slope_jump.py --linear --pmax 0.7 --no-axis) the cusps
  of n=3000 trace one smooth curve in D from ~1.9 near p*=1/2 down to ~0 at p*=0.652, EXCEPT a
  handful sitting near the floor at D ~ 0.03-0.08.  Investigated as possible computation errors.
- They are genuine.  Independent 50-digit check (exact binomials, no TINY, no recurrence, both
  rankings built directly): E'_- < 0 < E'_+ and D agree with the double result to 5 digits, and
  E(p*) is a local minimum by finite differences at steps 1e-10..1e-12.  (A first attempt with
  step 1e-7 wrongly said "not a minimum": at n=3000 the tie points are ~2e-7 apart, so that step
  crosses other tie points and measures nothing local.  Recorded so nobody repeats it.)
- Two populations among the 13 cusps of n=3000 with D < 0.2:
    11 are the F3<0 cusps: WIDE ties, j-i ~ 165-186, band i+j-n ~ 630-910, pair mass
       f(i) ~ 3.6e-5 .. 1.1e-4, F3 ~ -1700 .. -2000.  Every F3<0 cusp at n=3000 has D < 0.079,
       against a median D of 1.88 for F3>0 cusps (width median 43, f(i) median 1.1e-2).
     2 are the opposite: the NARROWEST possible ties, j-i = 1 and 2, at the very top of the cusp
       range (p* = 0.6519, 0.6521; max cusp p* is 0.65212).  Big mass (1.5e-2) but D is small
       because D = (j-i) f/(p*q*) and j-i is 1 or 2.
- So "F3<0 cusp" and "cusp with a tiny kink" are the same set, up to the two edge ties.  This
  matches section 3 (F3<0 cusps are shallow dips of wide ties) from the slope side: a wide tie has
  a small f(i) hence a small kick, and can only be a cusp where E' is already within ~D of zero --
  i.e. just past a smooth maximum, which is exactly the "always close to another cusp / narrow
  near-mode tie" observation.  Not analysed further.

### 2026-09-21 (Claude Code): the cusp mass floor is a rare-event statistic, not computation noise
Asked whether the wild scatter in the per-n minimum cusp mass could be numerical.  It is not.
- The 12 smallest-mass cusps in the whole n<=3000 table were recomputed at 60 digits with exact
  binomials, no TINY and no recurrence.  All are genuine cusps; f(i) agrees to ~1e-12 relative and
  S_- to ~1e-10 absolute.  Record: n=1075, (572,734), f = 6.562e-8.
- They ARE close to the decision threshold, which is worth knowing: a cusp of mass f needs
  |S_-| < kappa = (j-i) f, so f ~ 1e-7 forces |S_-| ~ 1e-5, within an order or two of MARGIN=1e-6.
  The record case has |S_-| = 1.571e-6, only 1.6x MARGIN -- but its absolute error is 3.3e-11, so
  the verdict still holds with ~5 digits to spare, and anything under MARGIN escalates anyway.
  Across the whole table 91 cusps have |S_-| < 1e-5 and only 44 needed interval arithmetic.
- WHY it scatters: kink_pos = S_-/kappa is UNIFORM on (-1,0) among cusps.  Decile counts at
  n=500/1000/2000/3000 give chi2 = 1.7, 1.2, 1.4, 2.3 on 9 dof (5% critical value 16.9) -- as flat
  as it is possible to look.  So "is this tie point a cusp?" behaves exactly like "does S_- land in
  a window of width kappa", and since kappa is proportional to f, low-mass tie points are cusps with
  probability proportional to their mass.  The per-n minimum is therefore a RECORD over ~n^2/4
  trials with a probability that vanishes linearly in f, which is precisely the kind of extreme-value
  statistic that scatters by orders of magnitude between neighbouring n.  Nothing to fix.
- A predictive model falls out, and it is the right shape but not yet calibrated:
  E[#cusps with mass < F] ~ sum over tie points with f<F of kappa/(spread of S_-).  With the IQR of
  S_- as the spread it over-predicts by ~3-6x (n=1000: 9.4 predicted vs 3 observed at F=1e-4;
  n=3000: 48.8 vs 8), because S_- is not uniformly spread over its IQR near zero -- the density of
  tie points near a zero of E' is what actually matters.  Getting that density right would turn
  the observed floor into a predicted distribution, and is the most promising route to a real
  bound: not "prove f > c" but "P(a cusp with f < F exists at this n) < epsilon".

### 2026-09-21 (Claude Code): cusps ARE the upward zero crossings of E'
Pursuing the density of S_- near zero, to calibrate the low-mass cusp model.  The density framing
turned out to be the wrong picture, and the right one is exact:
- E' DECREASES between tie points and jumps UP by D at each one (fact 11).  Therefore E' can only
  cross zero upward AT a tie point, and such a tie is by definition a cusp; it can only cross zero
  downward BETWEEN tie points, which is a smooth local maximum of E.  Cusps and smooth maxima
  strictly alternate.
- VERIFIED exactly: cusps == upward crossings, 181/181 at n=500, 352/352 at n=1000, 1061/1061 at
  n=3000; downward crossings number one fewer each time (E' starts negative just above p=1/2 and
  ends positive).  E'_-(t+1) <= E'_+(t) at all but 1 of 2,247,000 consecutive pairs at n=3000.
- So the cusp COUNT is structural -- the number of oscillations of E' about zero -- not a matter of
  chance.  That explains why the earlier "probability that S_- lands in a window of width kappa"
  model over-predicted by a steady ~2.9x: consecutive tie points can all satisfy E'_- in (-D,0),
  but only the crossing is a cusp.  The uniformity of kink_pos (previous entry) is a statement
  about WHERE in the jump the crossing falls, which is a different and still-open question.
- The crossing is size-biased toward large D: the first tie whose jump carries E' over zero wins, so
  cusps have median D = 1.88 while typical tie points have D astronomically smaller.  A LOW-mass
  cusp requires E'_- to have arrived within D of zero already.
- Modelling the per-n minimum D as Exponential (linear tail, independent crossings) FAILS: the
  probability-integral transform over n=200..3000 gives chi2 = 174 on 9 dof (5% critical 16.9),
  skewed so that observed minima are consistently LARGER than predicted.
- The reason is that the low-D "tail" is not a tail at all but a SECOND POPULATION, and it is
  exactly the F3<0 cusps.  Pooled over n<=3000: F3<0 cusps have D in [0.0000, 0.0816] with median
  0.0476; F3>0 cusps have D in [0.0267, 2.0913] with median 1.8828.  Every cusp with D<0.05 is F3<0
  (9118 of 9119), and the two groups barely overlap.  Fitted tail: count(D<D0) ~ D0^1.14 for the
  F3<0 group, not the D0^1.00 a uniform-arrival model predicts.
- So the mass floor is set by the F3<0 population specifically, and a usable probabilistic bound
  needs that population's own structure -- not a generic crossing argument.  That ties the floor
  question directly to the open F3-sign question rather than being independent of it.

### 2026-09-21 (Claude Code): kernel windowed, range widened to 0<=i<j<=n
Two changes to the certified kernel, done in this order so the first could be validated against
the existing data before the second changed what that data contains.
- WINDOW.  Only masses at or above TINY matter; the rest are zeroed, contribute exactly 0.0 to
  every sum, and sit as an equal block at the bottom of the ranking.  The masses fall away
  monotonically from the mode, so the recurrence stops once it drops below TINY and every pass runs
  over that window instead of [0,n].  Byte-identical by construction (adding exact 0.0 changes no
  running sum, Neumaier compensation included) and verified so: rebuilding n=3..1000 reproduced all
  998 per-n files and the interval-check log exactly.  Measured screening speedup 1.09x at n=1000,
  1.36x at 2000, 1.60x at 3000 -- real but below the window fraction (74%/49%/31% of masses),
  because the per-tie overhead does not shrink.  One trap handled: f[j]=f[i] is applied AFTER
  zeroing and can rescue a j just under TINY, so the recurrence refuses to stop before reaching i
  and j when f(i) is above TINY.  In 12,000 sampled ties i and j were never on opposite sides.
- RANGE.  j<=n now.  Counts go from n^2/4 to n^2/4 + (n-1).  Validation: 997 of 998 per-n files for
  n<=1000 byte-identical; the only change is n=3 gaining (1,3) at p*=0.6339745962, the one cusp this
  family contains anywhere.  Nine new interval checks, all j=n, all NOT.  cusps_all.csv is now
  1,591,533 cusps.
- analyze_cusps.tie_points() also regenerated tie points with the old range.  Fixed: it is used for
  nearest-neighbour gaps and intervening-tie counts, and the (i,n) points sit between the others in
  p* order, so excluding them would have corrupted every such metric.
- Parquet partitions rebuilt, n = 100..1000 by hundreds plus 2000..8000 by thousands.  Row counts
  rise by exactly n-1; cusp counts are unchanged (n=1000: 353, n=2000: 709, n=3000: 1062), while
  certified counts rise (n=3000: 247 -> 269) from the new j=n pairs needing interval arithmetic.

### 2026-09-22 (Claude Code): error analysis for the sharpened trigger -- MEASURED, and a warning
Step 1 of the sharpened CHECK trigger was to replace the asserted delta_f = 1e-11 and
delta_p = 1e-12 with derived bounds.  Derived the FORM from the algorithm, then measured every
quantity against 50-digit values (scale-free mpmath recurrence, no big binomials).

MEASUREMENTS (max over 20-30 random tie points per n; T denotes sum_k |w_k f_k (k - n p*)|):
- delta_p, relative error of the computed p*: <= 5.8e-12 (worst case n=8000, m=j-i=1).  The derived
  form is delta_p ~ q * eps * n * ln(n) / m, which at n=8000, m=1 gives 8e-12 against 5.8e-12
  measured -- form confirmed, and my first algebra was ~7x pessimistic.  Scales as 1/m, so wide ties
  are far better.  p* is NOT the error driver: its contribution to S_- is negligible.
- delta_f, relative error of a computed mass: ~1e-12 at n=2000, growing roughly as eps*n^1.4.
  So the asserted 1e-11 holds to about n=10,000 and fails beyond; any real bound must be
  n-dependent.  Error is largest at the window edge, smallest at the mode, as the recurrence
  accumulation predicts.
- error of S_- itself: 1.4e-10 (n=500), 8.5e-10 (1000), 3.0e-9 (2000), 1.2e-8 (3000), 3.5e-8 (5000).
  In units of eps*T this is 150, 315, 389, 866, 1144 -- i.e. observed err ~ 0.65 * eps * n^0.88 * T,
  and T ~ n^1.52.
- A first attempt at a bound, [(n+1)eps + 40 eps n^1.5] * T, over-covered the true error by
  950-4250x and would have flagged everything.  Recorded because the failure mode is the point:
  a bound calibrated by guesswork is useless even when it is valid.

THE RIGOROUS BOUND, and a correction to two claims I made while deriving it.
The bound on the error of S_- has two terms, T = sum|w f a|:
      summation (Sm is accumulated PLAIN; E is compensated)   (n+1)*eps*T
      mass error (each f carries relative error <= delta_f)   delta_f*T,  delta_f ~ 0.18*eps*n^1.4
Against MARGIN = 1e-6:
      n        T        summation    mass term    total     MARGIN/total
      1000  1.21e+04   2.69e-09     7.66e-09     1.04e-08      97
      2000  3.47e+04   1.54e-08     5.80e-08     7.34e-08      13.6
      3000  6.43e+04   4.28e-08     1.90e-07     2.32e-07       4.3
      5000  1.40e+05   1.55e-07     8.42e-07     9.97e-07       1.0
      8000  2.85e+05   5.07e-07     3.32e-06     3.83e-06       0.3
CORRECTION 1: I first reported MARGIN/bound using the summation term ALONE, giving 23 at n=3000 and
"exhausted near n=12000".  Including the mass term, which DOMINATES by ~6x, the true figures are 4.3
at n=3000 and exhaustion near n=5000.
CORRECTION 2: I proposed compensating the Sm accumulation as the fix, and it is not.  Implemented and
measured: it changes the OBSERVED error not at all (old 4.12e-09 total over 7 sampled ties, new
4.24e-09) because the mass error dominates the observed error by ~140x, and it reduces the RIGOROUS
bound only from 3.83e-06 to 3.32e-06 at n=8000 -- 13%, on a bound that exceeds MARGIN either way.
It costs 6-8% runtime and breaks byte-identity with every existing file.  REVERTED.
WHAT IS ACTUALLY TRUE: the rigorous bound is loose by ~2000x (3.3e-06 against an observed 1.7e-09 at
n=8000), because delta_f*T assumes every mass error aligns adversarially in sum_k w_k (df_k) a_k,
and they do not.  The data is also nowhere near the threshold: the smallest |S| ever decided in
double is 1.2e-05 at n=3000 and 1.9e-06 at n=8000, so the real safety factors are ~280 and ~1150
against the observed error.  A USEFUL rigorous bound needs the cancellation in that sum accounted
for -- that is the open piece of work, not compensation and not a bigger MARGIN.

CONSEQUENCE for the trigger: its win does NOT come from replacing MARGIN.  At n=1500 and n=2000
every single CHECK came from the near-tie proxy and none from the margin, and the margin term only
gets tighter than MARGIN below n~3000 while getting looser above n~8000.  The win comes from
replacing GAP = 1e-8 with the re-ranking bound, whose ambiguity threshold is ~delta_f ~ 1e-12 --
about 5000x tighter, which is what removed 99% of the checks in the n<=1000 prototype.  So the
trigger should be built around the cluster bound, and the margin term should be the rigorous one,
not a fitted constant.

### 2026-09-21 (Claude Code): the sharpened CHECK trigger -- IMPLEMENTED and VALIDATED
Built as the previous entry concluded: keep MARGIN, replace GAP.  binom_core gains _err_bounds()
and a second tag; both triggers are computed in the same pass, `sharp` picks which one drives the
verdict.  DEFAULT SINCE 2026-09-21: the sharpened trigger, in core.screen and in cusps_fast.py;
cusps_fast.py --legacy-trigger puts the GAP proxy back, and validate_trigger.py passes sharp
explicitly so it keeps comparing the two whatever the default is.  Re-validated after the flip:
same 195,062/181/0, same 0 newly flagged, and n=1400..1450 generated with the new default is
byte-identical to the same range generated with --legacy-trigger (62.4 s vs 73.2 s on 8 workers).

THE BOUND.  Swapping adjacent ranks of masses a,b moves S_- by f_a(a-np*) - f_b(b-np*).  So masses
whose order double precision cannot resolve are grouped into maximal clusters; within a cluster of
size c any rank moves by at most c-1, giving (c-1)*sum_{k in C} f_k|k - n p*|, summed over clusters.
That bound is added to MARGIN on BOTH sides (it perturbs S_- and S_+ equally; kappa is exact), so
no new constant enters.  Two masses are unresolvable when their relative gap is within
  2 df0 + dstep(|k-md| + |l-md|) + dp|k-l|/q
with df0 = 4 eps L + eps (L = |lnC[md]| + |md lnp| + |(n-md)lnq|; exp turns an ABSOLUTE argument
error into a relative one and the three n-sized terms do not cancel), dstep = 3 eps per recurrence
step, and dp = q*3 eps max(lnC[i],lnC[j])/m the relative error of p* (the lnC[i]-lnC[j] cancellation
costs a factor n/m).  Every constant counts flops; nothing is fitted.  Normalisation by the mass sum
is deliberately excluded -- it is common to all masses and cancels in every ratio tested.
The cluster {i,j} alone contributes NOTHING: f(i)=f(j) exactly and w_j = w_i - 1 is the left-limit
ranking by definition, not a numerical guess.  A third mass joining them makes it count in full.

VALIDATION (validate_trigger.py, both directions; the previous entry's two open caveats).
- Against every tie point the GAP rule escalated for n<=3000 -- all 195,243, verdicts already on
  record in cusps/interval_checks.log: 195,062 decided in double with the SAME verdict mpmath
  certified, 181 still CHECK, 0 disagreements.  Of the 181, 180 have rbnd = 0 and min|S| < MARGIN,
  i.e. they are honest margin cases; exactly ONE is a genuine re-ranking flag (n=2590, i=791,
  j=2243, rbnd = 0.44 against S_- = -7.5e-3, certified NOT) and it comes from a near-mode cluster,
  precisely the situation the feasibility entry predicted could arise.
- The untested direction, now tested: over EVERY tie point of n = 135, 400, 800, 1000, 1100, 1200,
  1500, 2000, 2500, 3000, 4000, 5000 (~14.6M tie points), the new bound flags NOTHING the old rule
  decided.  Checks per n: 12->0 (1000), 87->0 (2000), 269->0 (3000), 668->1 (4000), 1411->1 (5000).
  The bounds on the dropped ones are ~1e-287: the GAP proxy was firing on deep-tail masses just
  above TINY, exactly as diagnosed.
- End to end, n=1400..1450 generated both ways: all 51 per-n CSVs BYTE-IDENTICAL (certified_by
  included -- every escalation there resolved to NOT, so none reached the output), interval checks
  1306 -> 2, wall clock 71.1 s -> 60.5 s on 8 workers.

COST AND PAYOFF.  Computing both triggers costs 7.3% on the screening kernel (n=3000: 51.0 s ->
54.8 s single core); dropping the old one would recover little, and keeping it is what makes the
regression test possible.  Measured per-n cost at n=3000: 51.0 s screening + 247 checks * 143 ms =
86.3 s, against 55.2 s sharpened -- 1.56x.  EXTRAPOLATED from those measured pieces (screening
~n^2.6, checks ~n^3.0, certify ~n^1.05): 2.6x at n=5000, ~4x at n=8000, ~5x at n=10000.  The win
grows because certification was overtaking screening.

THE ONE GENUINE RE-RANKING FLAG, dissected (n=2590, i=791, j=2243, certified NOT).
Worth recording because it is the whole mechanism in one example, and because it is NOT a deep-tail
artefact.  p* = 0.5974588782672404, n p* = 1547.42.  The ambiguous cluster is k = 1530 and k = 1565,
two of the LARGEST masses in the distribution (f = 0.01250077796 each), straddling the mean at
-17.42 and +17.58.  Masses near-symmetric about n p* are near-equal, so this is structural, not
accidental: the deviations sum to only 0.163.  Their true relative gap is 1.31e-12 (60 digits);
double computes 2.84e-12 against a resolution threshold of 3.12e-12, i.e. 91% of it -- double gets
the order right but has no way to know that.  Swapping them moves S_- by exactly 0.4375, and here
the cluster bound EQUALS that exact effect (ratio 1.00) because the two sit on opposite sides of the
mode, so |a_k| + |a_l| = |a_k - a_l|.  Against S_- = -7.47e-3 the possible perturbation is 59x the
quantity whose sign is wanted, so the trigger cannot decide and escalates.
BUT THE VERDICT WAS NEVER AT RISK, and this is worth stating because it shows where the trigger is
conservative.  There are only TWO candidate rankings here, and both give NOT: S_- = -7.4696e-03 as
computed, S_- = +4.3006e-01 swapped.  Neither is a cusp, because the perturbation moves S_- and S_+
TOGETHER (kappa is exact) while kappa = 1452*f(791) = 4.9e-195 is effectively zero, so the cusp
window S_- < 0 < S_- + kappa has width 5e-195 and neither candidate lands in it.  The trigger bounds
the perturbation rather than enumerating the candidate rankings, so it cannot see that; a trigger
that evaluated the (here: two) reachable rankings would decide this in double.  That refinement is
real but was not built: one tie point in 195,243 does not pay for it.
GENERAL POINT: a rank ambiguity can only flip a verdict when kappa is comparable to the ambiguity.
kappa = (j-i) f(i) is tiny exactly when the tied pair is in the tail, which is when the cusp window
is narrow and the answer is NOT regardless.  The dangerous combination is a LARGE kappa (tied pair
near the mode) together with a large near-mode cluster, and no such case has been seen.
NOT A DOUBLE TIE.  The cluster is itself the tie point p*(2590,1530,1565) = 0.597458878267249396,
which sits 9.0e-15 (relative 1.5e-14) from p*(2590,791,2243) = 0.597458878267240398.  Two DISTINCT
tie points agreeing to 14 digits, not one p* shared by two pairs.  The standing convention (flag any
double tie found) is not violated, but this is the closest approach seen so far and is the mechanism
by which one would announce itself.

STATUS AND WHAT IS NOT DONE.  Zero disagreements is necessary, not sufficient; correctness rests on
_err_bounds, which is derived but assumes lgamma and log are <= 2 ulp and exp <= 1 ulp.  The cluster
bound is loose for c = 2 with both masses on the SAME side of the mode, where the exact change is
|f_k a_k - f_l a_l| rather than f_k|a_k| + f_l|a_l|; tightening that is available but bought only one
tie point in 195,243, so it was left alone.  Existing tables need no regeneration -- the verdicts are
identical, so n<=3000 stands as generated.

### 2026-09-22 (Claude Code): cusp tables to n<=5000, and the data policy rebuilt
THE RUN.  n=3001..5000 generated in 9 h 08 m (32,878 s) on 8 workers, resuming from n<=3000; the
projection from the measured screening curve (5.09e-8 * n^2.595 single core, 7.16x on 8 workers) was
9.2 h.  n<=5000 complete: 4,421,154 cusps, 47,299 with F3<0 (1.07%), 0 UNRESOLVED.  The n<=3000
portion is byte-identical to the archived table.  First range generated with the sharpened trigger
as the default: 399 interval checks over 2000 values of n, and only 94 rows in the entire 4.4M table
were certified by mpmath rather than double.  Cusps per n is almost exactly LINEAR, 0.354*n.
Sizes: cusps/ 672 MB of per-n CSVs + 672 MB merged; cusps_n5000.csv.gz 292 MB.

ANALYSIS (analysis/summary.txt; NOT yet interpreted -- see the flag in section 5):
- max cusp p* FLAT at 0.6521-0.6523 from n=1250 to 5000, overall max 0.65693.  No drift toward 0.66.
- F3<0 rate pinned at 1.07% in every n-range, including the 2.83M new cusps.
- The F3<0 signature is unchanged: median n*gap 0.045 vs 0.499 for F3>0, median intervening ties
  101 vs 1109, median percentile-at-n 27 vs 51.
- Closest cusp to E(1/2) still ALWAYS in the first band (i+j=n+1); min (E-E(1/2))*n fell to 0.0068.

THE DATA POLICY, rebuilt (user decisions 2026-09-22).  Three tiers became two KINDS, divided by
shape rather than size, because shape determines growth.  RESULT = one row per n, linear, may track
the frontier.  SAMPLE = one row per cusp, n^2, PINNED.  Full statement in CLAUDE.md.
- GitHub Releases DROPPED.  All 5 releases deleted (762 MB of assets).  Every asset across all of
  them had recorded ZERO downloads, and half were duplicates (v1.1 and v1.2 each carried full copies
  of the same three archives).  Tags kept as code snapshots.  The repo now offers the code and the
  knowledge plus a pinned sample; the frontier stays local.
- PINS SET AND FINAL: DECADE = 100,200,500,1000,2000,3000,4000,5000 (extending it to include 4000
  and 5000 cost a measured 215 KiB of clone size); F3_NMAX = 3000.
- negF3_neighbors.csv RENAMED to cusps_F3_negative.csv -- the old name read as "the list of all
  cusps" when it is the F3<0 subset, 17,081 rows out of 1,591,533 (1.07%).  Content unchanged.
- data/manifest.csv dropped from git: a local build log (wall-clock seconds, build timestamps)
  indexing Parquet dumps nobody outside can obtain.  Still written locally.
- Size bands (<=1 MB fine / 1-2 note it / 2-5 human discretion / >5 hard stop) are now ENFORCED by
  sys.exit in make_public_data.py, not just documented.  It refused cusps_F3_negative.csv on its
  first run and demanded --allow-large from a person, which is the intended behaviour.

A MEASUREMENT THAT CORRECTS THE OLD POLICY.  It claimed "each new version costs its full size in
history forever".  That is wrong by ~5x: git stores blobs zlib- AND delta-compressed, and these
tables GROW BY APPENDING, so consecutive versions share nearly everything.  Measured: 12 versions of
the three public tables, 10.65 MB of raw content, pack into 1.95 MB; the whole repo with full
history is 2.3 MB to clone, of which code and prose are 0.22 MB.  Committing the 7.67 MB n<=5000
F3 table would have added 1.96 MB, not 7.67 MB.  Never reason about history cost from raw file
size -- measure it with a probe clone and `git gc`.

LINE ENDINGS.  All four CSV writers (analyze_cusps.py x2, make_public_data.py, dump_ties.py) now
pass lineterminator="\n".  Python's csv module defaults to the 'excel' dialect, which is CRLF on
EVERY platform regardless of os.linesep -- that is RFC 4180 conformance, not a platform bug, and it
is deliberately NOT environment-dependent (output must be reproducible across machines, which is
what every byte-identity check in this project relies on).  All four committed CSVs had silently
been CRLF while the multi-GB generated tables were LF, because cusps_fast.py formats rows with
f-strings instead of using csv.  That mismatch made `cmp` between a public file and its source
table fail on 2407 lines that were numerically identical to the last digit.

### 2026-09-23 (Claude Code): NO DOUBLE TIES for n<=8000 -- exhaustive and EXACT
Question: can two different pairs (i,j) and (k,l) share a tie point p* at the same n?  Answer, over
all 42,674,665,999 tie points of every n from 3 to 8000: no.

WHY A FLOAT CHECK IS NOT ENOUGH, and this is the point of the entry.  The closest DISTINCT tie
points found are 2.2e-16 apart (n=7329) -- one ulp -- and 4.4e-16 at n=6000, 6.7e-16 at n=5331.
Meanwhile the computed p* itself carries error up to ~4e-12 for narrow pairs, because
lnrho = (lnC[i]-lnC[j])/(j-i) cancels two numbers of size ~n ln2 and the cancellation costs a factor
n/m.  So beyond n~2000, "the gap is zero" and "the gap is tiny" are indistinguishable in double, and
a pure float scan would be reporting a confidence it does not have.

THE EXACT TEST.  rho^(j-i) = C(n,i)/C(n,j), so (i,j) and (k,l) collide iff
    (C(n,i)/C(n,j))^(l-k) == (C(n,k)/C(n,l))^(j-i),
an identity between integers.  Do NOT form them (C(5000,2500) has ~1500 digits, raised to a power of
thousands).  Take p-adic valuations: it holds iff for EVERY prime p<=n
    (l-k)(v_p C(n,i) - v_p C(n,j)) == (j-i)(v_p C(n,k) - v_p C(n,l)),
with v_p C(n,k) = (S_p(k) + S_p(n-k) - S_p(n))/(p-1) by Legendre, S_p = digit sum in base p.  No big
integers, and almost every candidate dies on the first prime.  The first 12 primes are tested
vectorised in numpy over all candidates at once; only survivors reach the scalar sweep.

WHAT WAS CHECKED.
- The 17 complete tie-point dumps (n=100..8000), 51.7M tie points: 1.35M candidate pairs below the
  1e-9 screen, all proved distinct.
- Cusp-cusp, every n<=5000, 4.4M cusps: minimum gap 1.489e-09, i.e. not one pair was even a
  candidate.  Cusps are far better separated than tie points in general.
- EXHAUSTIVE: every tie point of every n<=8000, 4.27e10 of them, 511 s on 8 workers.  Candidates
  below the screen grow as ~n^4 (6 at n=501, 2415 at n=2001, 626,191 at n=7995; ~1e9 in total).
  Zero collisions.  An independent single-threaded run to n<=5000 agreed.
The screen at 1e-9 is ~250x the worst-case numeric error in p*, so it cannot miss a true collision.

TOWARD A PROOF.  The valuation identity is probably the right entry point: it turns "two tie points
coincide" into "two distinct pairs have proportional binomial-coefficient exponent vectors", a
statement about factorisations rather than about real numbers.  Not attempted.

PERFORMANCE NOTE.  The first exhaustive run was single-threaded and the user observed it sitting on
the efficiency cores with the performance cores idle.  Jobs launched by the agent inherit a reduced
QoS; taskpolicy -c only clamps downward so it cannot be promoted from inside.  A worker pool sidesteps
it -- the same n<=1500 scan is 15.7 s serial and 3.0 s on 8 workers (711% CPU).  Always parallelise
long compute rather than trying to raise its priority.

### 2026-09-23 (Claude Code): how to SEARCH for tie-point collisions -- a prime obstruction
The exhaustive check (previous entry) settles n<=8000 but costs ~N^3 cumulatively, so it stops
there: N=50,000 would be 35 h and N=100,000 about 11 days, and the float pass also hits a memory
wall near n~25,000 (it builds arrays of n^2/4 entries -- 15 GB at n=50,000).  This entry is about
searching much further.

WHAT RESTS ON WHAT -- read this before any of the timings below.  There are three tiers, and they
are NOT three levels of confidence in the same method; two of them share one assumption exactly.
  TIER 1, proof-grade.  Examine EVERY pair with the exact simplifying test, then apply Facts A and
    B.  Both are proved, so a clean run proves no collisions in the range scanned, with no further
    assumption.  This is what the n<=8000 result rests on.
  TIER 2, rests on Fact C and nothing else.  Examine only pairs with width <= g(j) or band <= g(j),
    where g(j) = j - prevprime(j).  Exhaustive IF AND ONLY IF Fact C holds.
  TIER 3, extrapolation.  Use a constant cutoff (e.g. min(width,band) <= 3, the largest value ever
    observed).  Fast, justified by nothing, and it can silently miss a collision.
A BLANKET cutoff C = max prime gap below N is NOT a fourth tier: it is Tier 2 applied loosely.  It
carries exactly the same dependence on Fact C and is ~9x slower than the per-j form, because the
mean gap is 10.0 where the maximum is 113 (j<=10^6).  Do not present it as the safer option -- it
buys no certainty, only cost.  It is recorded here only so nobody reinvents it.

      COST TO SCAN EVERY n <= N, 8 workers
      N          TIER 1 (proved)   TIER 2 (needs Fact C)
      10,000     36 min            15 s
      100,000    25 days           34 min
      1,000,000  25,000 days       3.0 days
So Fact C is worth about a factor of 1000 at N = 100,000.  That is the whole reason to write it out
carefully rather than lean on it informally.

NOTATION.  For a tie point (i,j) write the ODDS VALUE r = p*/(1-p*), so p* = r/(1+r) and two tie
points collide exactly when their r agree.  Write WIDTH = j-i and BAND = i+j-n.  The defining
equation is  r^width = C(n,i)/C(n,j).  So r is a width-th root of a fraction.

FACT A (PROVED, algebra).  For a FIXED width, r is strictly increasing in i, so two tie points of
the same width can never collide.  Proof: C(n,i)/C(n,j) = prod_{t=i+1}^{j} t/(n+1-t), and
t/(n+1-t) is strictly increasing in t; sliding the window one step right drops the smallest factor
and adds a larger one, so the product strictly increases.  Also checked numerically at several n.

THE ROOT MAY SIMPLIFY.  x^6 = 4 is really x^3 = 2: a sixth root wearing a disguise.  Define the
TRUE ROOT ORDER M = the least M>0 with r^M rational.  M always divides the width.  Concretely, with
e_p the exponent of prime p in C(n,i)/C(n,j) and g = gcd_p e_p, we have M = width/gcd(width, g).
Call a tie point SIMPLIFYING when M < width.  Measured: 314 of the 5,353,299 tie points with
n<=400 simplify, about 1 in 17,000.

FACT B (PROVED, a deduction from Fact A).  M is a property of the NUMBER r, so colliding tie points
share it.  If neither simplifies then M = width for both, hence equal widths, which Fact A forbids.
THEREFORE ANY COLLISION HAS AT LEAST ONE SIMPLIFYING PARTNER.  That is the reduction: stop comparing
pairs, and instead enumerate the simplifying ones and look for each one's partner.  The partner
search is cheap -- it must have the same M and the same reduced fraction, and if it does not itself
simplify its width IS M, so Fact A makes it a binary search over that width.

FACT C -- A PRIME OBSTRUCTION.  WRITTEN OUT AND TESTED 2026-09-23 (see the entry of that date for
the full statement, the four tests, and the screen built on it).  It holds.
  C(n,i)/C(n,j) = j! (n-j)! / ( i! (n-i)! ).
  Let P be a prime with  max(i, n-i, j/2) < P <= j.  Then v_P(j!) = 1 (because j/2 < P <= j),
  v_P(i!) = 0 (P > i), v_P((n-i)!) = 0 (P > n-i), and v_P((n-j)!) = 0 (n-j < n-i < P).
  So v_P of the ratio is EXACTLY 1, hence gcd_p e_p = 1, hence M = width: NOT simplifying.
  The interval (max(i, n-i, j/2), j] has length  min(width, band, ceil(j/2)).
  CONCLUSION: a simplifying tie point requires that interval to be PRIME-FREE, so
      min(width, band)  <=  g(j) := j - (largest prime <= j).
  This EXPLAINS the observed data rather than fitting it: over n<=400 the 314 simplifying points
  have min(width,band) = 1 (286 of them), 2 (26) or 3 (2), and never more -- prime-free runs of
  length 1-3 are common below 400, longer ones are not.

WHY THIS MATTERS, AND A CORRECTION.  I first proposed identifying "the family" the 24 non-band-1
simplifying points belong to.  The user pointed out that this is worthless: a family fitted to
n<=400 says nothing about a new family at n=60,000, so the screen would stay an extrapolation with
a better name.  That objection is correct.  Fact C is different in kind -- it is a reason a
simplifying point CANNOT exist outside the prime gaps, not a description of where they have been
seen.  If the sketch holds up, the screen becomes exhaustive rather than heuristic.

THE SCREEN, AND ITS COST.  Per n: sieve the primes to n, compute g(j) = j - prevprime(j) for each j,
and examine only pairs with width <= g(j) or band <= g(j).  NOTHING IS DOWNLOADED -- g(j) comes from
a local sieve, and the per-j gap is far tighter than any blanket maximum (mean g is 10.0 for
j<=10^6 against a maximum of 113).  Published maximal-gap tables would only be needed if we wanted
a blanket bound beyond sieving range, which this formulation avoids.
      N          pairs examined          8 workers      full n^2 scan for comparison
      10,000     5.7e8                   15 s           36 min
      100,000    7.8e10                  34 min         25 days
      1,000,000  1.0e13                  3.0 days       25,000 days
Measured throughput for the exact width-organised scan is 4.8M pairs/core-s (against 10.4M for the
existing float-screen check, i.e. the exact scan is ~2x slower per pair but uses O(n) memory instead
of O(n^2), so it does not hit the wall the float method does).

SIDE PRODUCT: A TABLE OF SIMPLIFYING TIE POINTS.  The screen's primary output IS that table; finding
collisions is a second step on top of it.  For n<=400 it has 314 rows, each carrying width, band,
the true root order M and the reduced fraction.  15 of them have M=1, i.e. p* is EXACTLY RATIONAL --
all of width 2, where the ratio is a perfect square: n=50 (48,50) gives r=35 and p*=35/36; n=289
(287,289) gives r=204 and p*=204/205; n=244 (241,243) gives r=99 and p*=99/100.  The other 299 have
M from 3 to 195, e.g. n=9 (2,8) is the cube root of 2 and n=259 (4,256) is the 42nd root of 2.
(Width-1 tie points are rational too -- p* = (i+1)/(n+1) -- but trivially so, with M = width = 1, so
they are not counted as simplifying.)

STATUS.
  PROVED:    Fact A; Fact B given Fact A; the band-1 identity C(n,i)/C(n,n+1-i) = (n+1-i)/i.
  SKETCHED:  Fact C -- the argument is short and the five conditions on P look right, but it has not
             been written out rigorously.  Verified against 314 cases, 0 counterexamples.
  MEASURED:  the 1-in-17,000 rate, the min(width,band) <= 3 observation, all the timings.
  OPEN:      whether Fact C survives a careful write-up.  If it does, this screen is exhaustive and
             n<=100,000 is half an hour.  If it does not, the screen is still a good search and a
             bad proof, and should be described that way.

### 2026-09-23 (Claude Code): the cusp mass floor over n<=5000 -- a new record, 18x lower
From plotting/mass_floor_linear.py (linear n axis, one pixel per n) over the complete table n=3..5000.
- New smallest pair mass at any cusp: n=3076, (1476,1781), p*=0.52951, f(i) = 3.630e-9, against
  6.56e-8 (n=1075) for n<=3000.  It is certified by interval arithmetic (iv50; S_- = -1.19e-7, below
  MARGIN), and --recheck at 50 digits agrees: S_- = -1.1947e-7, S_+ = 9.878e-7, cusp True.  It also
  has F3 < 0 (-1342.4).  Next smallest above n=3000: 4.43e-9 (n=4712), 5.15e-9 (n=4314); then a gap
  to 2.24e-8 (n=4227).  The running minimum is flat at 3.63e-9 from n=3076 to 5000.
- Refitted over n=300..5000: median ~ 0.592 n^-0.500 (log-log corr -1.000), unchanged.  Low tail
  still follows no law: minimum ~ n^-1.57 (corr -0.60), 1st percentile ~ n^-0.98 (corr -0.31).
  So the rare-event reading of the 2026-09-21 entry stands; the record simply moved.
- n=2 has no cusp: its only tie point, (1,2) at p*=2/3, has S_- = 0 exactly (E = 4p-3p^2 just below,
  2p just above).  n=3 is the first n with a cusp, so plots over n start there.

### 2026-09-23 (Claude Code): primality of n or n+1 does not show in the cusp statistics
prime_compare.py, n=100..5000.  Per n: cusp count, median and minimum of f(i).  Compared WITHIN a
parity, since every prime above 2 is odd and n+1 prime forces n even: n prime vs odd composite n,
and n+1 prime vs even n with n+1 composite.  Each statistic is detrended by the median of log s over
same-parity neighbours within +-100 (self excluded); two-sided permutation p over 20,000 shuffles.
      statistic   n prime vs odd comp.     n+1 prime vs even      residual SD
      count       -0.02%   p=0.72          -0.00%   p=0.99         1.0%
      median f    +0.06%   p=0.031         +0.00%   p=0.93         0.64%
      min f       +1.5%    p=0.84          +14.5%   p=0.037        152%  (log units)
- No effect of any practical size.  The two p~0.03 values are 2 of 6 tests (Bonferroni 0.008) and
  are small against the spread: 0.06% on a 0.64% SD, 14% on a 152% SD.  The min-f one weakens to
  p=0.15 for n>=1000.  Window +-30 gives the same picture.
- Window +-300 invents a count/median "effect" (p~0.01) for BOTH primality classes at once.  That is
  an artefact: a wide window mis-fits the trend's curvature at small n, where primes are densest.
  Keep the window narrow for any test of this kind.
- Parity itself (odd vs even, the control) also shows nothing in these three statistics, unlike the
  E - E(1/2) of the lowest cusp, which splits cleanly by parity.
- A prime-n strip under the linear mass-floor plot was tried first; the user found it not useful
  and it was dropped.

### 2026-09-23 (Claude Code): odd vs even n -- only the lowest cusp's height splits by parity
parity_compare.py, n=100..5000, over every per-n column of public/per_n_summary.csv plus the
median and minimum of f(i).  Sawtooth test: z(n) = s(n) - (s(n-1)+s(n+1))/2 cancels the trend, and
delta = (mean z_odd - mean z_even)/2 is the odd-minus-even offset; permutation p over 20,000.
- min_E_minus_Ehalf: odd n sit 41.24% HIGHER than even n, delta/scatter = 1.00, p < 5e-5.  This
  is the known split (lowest_cusp.py: E-E(1/2) times n^1.5 is 0.6828 odd vs 0.4834 even; ratio
  1.4125), recovered independently -- so the test does detect a real parity effect.
- Everything else is null, |delta| <= 0.025 of the scatter, all p >= 0.10: n_cusps (-0.007%),
  n_negF3, n_double_cusps, max_pstar, min_pstar, min_F3, median_slope_jump, min_slope_jump,
  median f (+0.001%), min f.  n>=1000 gives the same picture.
- So parity acts on how HIGH the lowest cusp sits above E(1/2), not on how many cusps there are or how their
  masses and slope jumps are distributed.

### 2026-09-23 (Claude Code): no double ties for n<=10000 (extends the n<=8000 result)
Tier 1 throughout -- every pair examined, every pair double precision could not separate decided
exactly by p-adic valuations.  Rests only on Facts A and B, both proved; no Fact C, no cutoff.
  83,345,832,499 tie points, every n from 3 to 10000, 1018 s on 8 workers.  ZERO collisions.
  Predicted 17 min from the measured 511 s at n<=8000 scaled by (10/8)^3; actual 1018 s.
  Roughly 3e9 candidate pairs fell under the 1e-9 float screen and were decided exactly.
FLOATS ARE NOW DEMONSTRABLY USELESS HERE, not just marginal.  At n=7500 and n=8333 the minimum gap
between distinct tie points is EXACTLY 0.0 in double precision -- two different (i,j) pairs produce
bit-identical p* -- and the exact check proves them distinct.  Any screen that stopped at the float
comparison would now be reporting collisions that do not exist.
COST OF GOING FURTHER at Tier 1 (~N^3): 2.4 h to 20,000, 8 h to 30,000, ~2 days to 55,000.  Memory
binds first: each worker holds n^2/4 entries, 3.75 GB per worker at n=25,000, so 8 workers need
~40 GB there (69 GB available).  Fewer workers trades speed for headroom.  The prime-gap argument
in the previous entry would replace this with ~34 min to 100,000, but only if it holds up.

### 2026-09-23 (Claude Code): FACT C written out, and the accelerated screen built on it
STATEMENT.  Tie point (i,j), 0 <= i < j <= n, i+j > n; width m = j-i, band b = i+j-n.  Since
C(n,k) = n!/(k!(n-k)!),
        R = C(n,i)/C(n,j) = [ j! (n-j)! ] / [ i! (n-i)! ].
  CLAIM.  Let P be a prime with  (1) P <= j,  (2) P > j/2,  (3) P > i,  (4) P > n-i.
          Then v_P(R) = 1 EXACTLY.
  PROOF.  v_P(j!) = floor(j/P) + floor(j/P^2) + ...  By (1)+(2), 1 <= j/P < 2 so floor(j/P) = 1;
          by (2), P^2 > j^2/4 >= j for j >= 4, so every higher term vanishes: v_P(j!) = 1.
          By (3) no multiple of P is <= i, so v_P(i!) = 0.  By (4), v_P((n-i)!) = 0.
          j > i gives n-j < n-i < P, so v_P((n-j)!) = 0.  Hence v_P(R) = 1+0-0-0 = 1.   []
  COROLLARY.  Such a P forces the point NOT to simplify: an exponent equal to 1 makes the gcd of
  all exponents 1, so d = gcd(m,1) = 1 and the true root order M equals the width.
  THE INTERVAL.  (1)-(4) say max(i, n-i, j/2) < P <= j, and max(i, n-i) < j always (i < j by
  definition; n-i < j is exactly i+j > n).
  CONTRAPOSITIVE.  If the point simplifies there is no prime in that interval.  Let P* be the
  largest prime <= j.  Bertrand gives P* > j/2, so the j/2 term never binds and we need
  P* <= max(i, n-i).  Therefore
        min(width, band) = j - max(i, n-i)  <=  j - P*  =:  g(j).
  SPECIAL CASE: j prime => P* = j => g(j) = 0 => the point can NEVER simplify.
  EDGE CASE: v_P(j!) = 1 needs j >= 4; j <= 3 is checked directly.

TESTS (all pass).
  1. No simplifying point among the 314 with n<=400 has j prime.                       0 violations
  2. All 314 satisfy min(width,band) <= g(j).                                          0 violations
  3. v_P(R) == 1 for a qualifying P, on 2460 random non-simplifying pairs.             0 wrong
  4. Bertrand numerically: prevprime(j) > j/2 for j = 4..2000.                         0 violations
  5. THE ONE THAT MATTERS -- screen_collisions.py --verify runs the RESTRICTED window and the
     UNRESTRICTED scan at every n and compares.  Identical at every n up to 1500 (1461 simplifying
     points).  So the window provably-and-empirically loses nothing in that range.

screen_collisions.py.  TIER 2: exhaustive iff Fact C holds.  Per n it examines only pairs with
width <= g(j) or band <= g(j) -- 39,246 candidates at n=4000 against n^2/4 = 4,000,000, a 102x
reduction.  Then Fact B: every collision has a simplifying partner, so for each simplifying point
it searches for the partner, by binary search over width-M pairs (valid by Fact A's monotonicity)
followed by exact exponent-vector comparison, plus a direct pairwise check among the simplifying
points themselves.  The binary search is never exercised by real data (no collisions exist), so it
was validated by PLANTING targets: 299 of 299 located.

PERFORMANCE, and three mistakes worth recording.  The first version was 207 s to n<=2000 and scaled
as n^3.2 -- years at 100,000.  Fixes, in order of payoff:
  (a) the candidate list was built with a Python loop over j; np.repeat made it ~50x faster;
  (b) the d>1 prefilter used a gcd over 8 primes, leaving ~2^-8 of candidates alive -- still ~80
      per n, each costing a full Python sweep over pi(n) primes, which was 94% of the runtime.  An
      ESCALATING filter (8 primes, then 32, then 160) cut that to almost nothing.  Valid because
      the gcd over any SUBSET of primes is a multiple of the true gcd, so gcd(m, G_subset) = 1 is
      always a sound rejection;
  (c) the later stages still built whole length-n valuation arrays for a handful of survivors;
      evaluating Legendre pointwise on the survivors only gave another ~1.5x.
Now: 4.9 s to n<=4000, 20 s to n<=8000, 59 s to n<=12000.
