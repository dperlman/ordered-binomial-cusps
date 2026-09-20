# Research log: the ordered binomial distribution and its cusp points

Written 2026-09-18 from a claude.ai project conversation (Sept 9–18, 2026). Everything below was
established in that conversation; "verified" means checked numerically with the stated range,
"proved" means an argument was given and checked. Keep appending to this file.

## 1. Definitions and notation

- f_p(k) = C(n,k) p^k (1-p)^(n-k), k = 0..n.  q = 1-p.  K ~ Bin(n,p).
- Ranks: w_k = rank of f_p(k) among the n+1 masses in increasing order, 0 = smallest.
- Ordered-binomial expectation:  E(n,p) = sum_k w_k f_p(k).
- G(p) = sum_{0<=k<l<=n} |f_p(k) - f_p(l)|.   Then  G = 2E - n.
- Tie point of the pair (i,j), 0<i<j<n:  p*(n,i,j) = the p with f_p(i)=f_p(j):
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
- Extend cusp tables to n<=5000 (running); check whether negative-F3 cusps remain "close" and
  whether max cusp p* stays ~0.65.
- Quantify cusp "depth" (dip height before the nearest smooth max) for all cusps; conjecture:
  negative-F3 cusps are the shallowest.
- Prove piecewise concavity of E (fact 7 above is numerical).
- Prove that the first tie point above 1/2 is the innermost pair (fact 5 is numerical for n<=200).

## 6. Files

- cusps_fast.py: the generator (numba); its docstring has the definitions and method.  The older numpy
  script cusps_parallel.py was removed 2026-09-18 after cusps_fast.py reproduced all of its n<=1000
  output byte for byte (see section 7).
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

## 7. Log entries

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
