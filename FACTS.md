# What is established

A quick list of what can be said for certain so far.  No proofs here -- each item points to where
the proof or the run is recorded in `RESEARCH_LOG.md`.  Notation as in `CLAUDE.md`: tie point
p\*(n,i,j) with 0 ≤ i < j ≤ n and i+j > n; width = j−i; band = i+j−n; E(n,p) the ordered-binomial
expectation; a *cusp* is a tie point that is a local minimum of E.

**Proven** = a mathematical proof exists.  **Screened** = established by an exhaustive or certified
computation over the stated range, by the named program.  Dates are when the result was locked in.

## Proven mathematically

| # | Statement | Locked in | Where |
|---|---|---|---|
| P1 | E(n,p) = n/2 + ½ Σ<sub>k<l</sub> \|f<sub>k</sub> − f<sub>l</sub>\|.  E(n,p) = E(n,1−p); E(n,0) = E(n,1) = n. | ≤ 2026-09-20 | Log §2, fact 3 |
| P2 | At a tie point, i+1 ≤ (n+1)p\* ≤ j, so i < np\* < j.  Equality only for width 1, where p\* = (i+1)/(n+1) exactly. | ≤ 2026-09-20 | Log §2, fact 1 |
| P3 | At a tie point the rank of the pair is w<sub>i</sub> = i + n − j. | ≤ 2026-09-20 | Log §2, fact 2 |
| P4 | Cusp test: a tie point is a cusp iff −(j−i)f(i) < S₋ < 0 (equivalently S₋ < 0 < S₊). | ≤ 2026-09-20 | Log §2, fact 6 |
| P5 | Every local minimum of E(n,·) is at a tie point: E has no smooth local minima. | 2026-09-21 | Log §2, fact 10 (external paper, via the user) |
| P6 | No tie point is a local maximum: E′ jumps **up** by D = (j−i)f(i)/(p\*q\*) > 0 at every tie point.  With P5: local minima are exactly the cusps; local maxima are smooth points. | 2026-09-21 | Log §2, fact 11 |
| P7 | Above the last tie point, (n−1, n) at p\* = n/(n+1), all masses are in natural order and E = np exactly. | 2026-09-21 | Log 2026-09-21, "the 0<i<j<n convention" |
| P8 | The inequality 2np\* < i+j is false in general (smallest counterexample n=9, (2,8)). | ≤ 2026-09-20 | Log §2, fact 9 |
| P9 | **Fact A.** At fixed width, p\* is strictly increasing in i.  So two tie points of the same width never coincide. | 2026-09-23 | Log 2026-09-23, collision search |
| P10 | **Fact B.** If two tie points coincide, at least one of them is *reducing* (its odds value p\*/(1−p\*) is a lower-order root than its width suggests). | 2026-09-23 | Log 2026-09-23, collision search |
| P11 | **Fact C.** A reducing tie point satisfies min(width, band) ≤ j − (largest prime ≤ j).  In particular it can never have j prime.  Independently checked by a second model 2026-09-24. | 2026-09-23 | Log 2026-09-23, "FACT C written out" |
| P12 | Band 1: C(n,i)/C(n,n+1−i) = (n+1−i)/i, so every band-1 tie point has (p\*/q\*)<sup>width</sup> = j/i, and p\* → ½ + 1/(2(n+1)). | 2026-09-23 | Log 2026-09-23 and 2026-09-24 |
| P13 | Counting.  There are n(n+1)/2 tie points (pairs 0 ≤ i < j ≤ n) in (0,1): ⌊n²/4⌋ above ½, ⌊n²/4⌋ below ½ (mirror images), and ⌈n/2⌉ exactly at p = ½.  A pair has p\* = ½ if and only if i + j = n, so all ⌈n/2⌉ of those coincide there.  Cumulatively, Σ<sub>n=3..N</sub> ⌊n²/4⌋ = ⌊N(N+2)(2N−1)/24⌋ − 1.  The number of **distinct** p\* values in (0,1) is ⌊n²/2⌋ + 1 exactly when there are no other collisions — so for n ≤ 100,000 (S2). | 2026-09-24 | Log 2026-09-24, "tie-point counts" |

## Screened computationally

| # | Statement | Range | Program | Locked in |
|---|---|---|---|---|
| S1 | **No tie-point collisions** — no two different pairs (i,j), (k,l) at the same n have the same p\*.  Confirmed by three independent methods. | 3 ≤ n ≤ 10,000 | `check_collisions.py --exhaustive`; `screen_collisions.py --no-fact-c`; `screen_collisions.py` | 2026-09-24 |
| S2 | **No tie-point collisions.** Relies on P9–P11. | 3 ≤ n ≤ 100,000 | `screen_collisions.py` | 2026-09-24 |
| S3 | Complete catalogue of reducing tie points: 163,529 of them (12,667 for n ≤ 10,000, where it is also confirmed without P11).  217 have p\* exactly rational beyond the trivial width-1 case.  None has min(width, band) > 4. | n ≤ 100,000 | `screen_collisions.py --save` (local file `simplifying_n100000.csv`) | 2026-09-24 |
| S4 | Complete certified cusp tables: 4,421,154 cusps, 0 unresolved.  Each decision certified (double-precision screen with margin, interval arithmetic where needed). | 3 ≤ n ≤ 5000 | `cusps_fast.py` | 2026-09-22 |
| S5 | Every cusp has p\* < 0.657 (maximum 0.65693, at n = 15).  For every n > 1250 the per-n maximum lies between 0.6517 and 0.6525. | n ≤ 5000 | `cusps_fast.py`, `analyze_cusps.py` | 2026-09-22 |
| S6 | E(p\*) > E(1/2) at every cusp — the global-minimum conjecture holds at every local minimum found.  (E in double precision; smallest margin ~1.4×10⁻⁶ against an error ~10⁻¹².) | n ≤ 5000 | `analyze_cusps.py` | 2026-09-22 |
| S7 | The lowest cusp is always in band 1.  It is also the first cusp (smallest p\*) except at n = 6, 21, 50, 76, 125. | n ≤ 5000 | `analyze_cusps.py`; scan of `cusps_all.csv` | 2026-09-23 |
| S8 | Band 1 contains exactly one cusp, except exactly two at n = 6, 9, 21, 50, 76, 125, 321. | n ≤ 5000 | scan of `cusps_all.csv` | 2026-09-23 |
| S9 | Pairs (i, n) hold no cusps for n ≥ 4; the only one is n = 3, pair (1,3). | 3 ≤ n ≤ 5000 | `cusps_fast.py` | 2026-09-22 |
| S10 | F3 < 0 at 47,299 cusps (1.07%), a rate flat in n. | n ≤ 5000 | `analyze_cusps.py` | 2026-09-22 |
| S11 | Cusps are exactly the upward zero crossings of E′; cusps and smooth maxima alternate. | n = 500, 1000, 3000 | log 2026-09-21 | 2026-09-21 |
| S12 | E is concave between consecutive tie points. | n ≤ 100 | pre-repository work | ≤ 2026-09-20 |
| S13 | The first tie point above ½ is the innermost pair. | n ≤ 200 | pre-repository work | ≤ 2026-09-20 |
| S14 | Number of cusps above ½ ≈ 0.3537 n (within about ±5 of the fit at any n, worst ±20); over all of (0,1), ≈ 0.7074 n + 1 (mirror images plus the cusp at p = ½).  Empirical fit, not a formula. | 500 ≤ n ≤ 5000 | fit to `public/per_n_summary.csv` | 2026-09-24 |
| S15 | Fraction of tie points that are cusps ≈ 1.415 / n (same leading term for the ⌊n²/4⌋ tie points above ½, all tie points in (0,1), or distinct p\* values).  Follows from S14. | 500 ≤ n ≤ 5000 | as S14 | 2026-09-24 |
| S16 | **(★) at every switch point**: D(p<sub>m</sub>) < D(½) for all n−1 switch points p<sub>m</sub> in (½,1).  Scaled margin 2(D(½)−D(p<sub>m</sub>))n<sup>1.5</sup> is smallest at m = 1 for every n and strictly increasing in m; minimum 0.2573 (n = 4), ≥ 0.2974 for n ≥ 100.  Double precision; worst measured error (mpmath, 50 digits) 4.6×10⁴ below the smallest margin.  Implies the main conjecture for these n **only together with the switch-point lemma**, which is argued, not proved. | 3 ≤ n ≤ 5000 | `star_check.py` | 2026-09-24 |

Screened results rest on the code being correct (S1–S3: externally reviewed and cross-checked
against independent brute force) and on standard library functions (`lgamma`) behaving as documented.

## Not established — open

- The main conjecture: E(n,p) ≥ E(n,½) for all p, all n.  (S6 is numerical support only.)
- The switch-point lemma (D quasi-convex between switch points; Log 2026-09-24, claude.ai entry, item 2): argued, not written up.  S16 needs it.
- That no cusp lies above p ≈ 0.66 for all n.  (S5 is numerical.)
- The sign of F3 at cusps; whether F3 < 0 cusps have a structural characterisation.
- A lower bound on the pair mass f(i) at a cusp.
- No tie-point collisions for **all** n.  (S1–S2 are bounded ranges.)
- Concavity of E between tie points in general (S12).
