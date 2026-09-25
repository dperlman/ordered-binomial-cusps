# Cusp points of the ordered binomial distribution

Tools and certified data for the **ordered binomial expectation**

$$
E(n,p) = \sum_{k=0}^{n} w_k f_p(k)
\qquad \text{where} \qquad
f_p(k) = \binom{n}{k} p^k (1-p)^{n-k}
$$

and `w_k` is the **rank** of the mass `f_p(k)` among all `n+1` masses, smallest = 0. Equivalently,
and without any sorting:

$$
E(n,p) = \sum_{k \lt l} \max(f_k, f_l)
       = \frac{n}{2} + \frac{1}{2} \sum_{k \lt l} \lvert f_k - f_l \rvert
$$

Because the ranking changes as `p` moves, `E(n,.)` is a piecewise-smooth curve: concave on each
piece, with a convex **kink** wherever two masses swap order. Those crossings are the **tie points**
`p*(n,i,j)`, where `f(i) = f(j)`, and the ones where `E` turns around — where the derivative changes
sign — are the **cusp points**. Every local minimum of `E` is one of them.

There are a lot of tie points (about `n²/4`) and rather few cusps (about `0.354n`, so roughly one in
700 at `n = 1000`), and deciding which is which is numerically delicate: the test is the sign of two
quantities that can sit `1e-7` from zero while the masses involved span 300 orders of magnitude.
**This repository is the machinery for getting those decisions right, at scale.**

## Two repositories

| | |
|---|---|
| **[OBDExplorer](https://github.com/dperlman/OBDExplorer)** | Interactive exploration — 2D/3D/4D plots, live views of the structure. Go there to *look at* the distribution and develop intuition. |
| **this repo** | Optimised, certified tools for *producing the raw data* those explorations rest on, plus the datasets themselves. |

## What is certified, and what that means

Every cusp decision is proved, not estimated. A fast double-precision screen decides the easy cases;
anything near a decision boundary, or where two masses are too close for double precision to order
reliably, escalates to **mpmath interval arithmetic** at 50, 100 or 200 digits. A verdict is only
recorded when the intervals separate. Across `n <= 3000` that path was taken 195,243 times, and
**no tie point anywhere is left unresolved**. As an independent check, `n <= 200` reproduces a
separate 50-digit computation exactly.

## Current results

- Certified cusp tables complete for `n <= 5000`: **4,421,154 cusps**, 47,299 of them with `F3 < 0`.
- Tie-point datasets (*every* tie point, not just cusps) for `n` = 100…1000 by hundreds, then
  1000…8000 by thousands — up to 16 million rows each.
- `E(p*) > E(1/2)` at **every one** of those 4.4 million cusps, which is the conjecture below
  holding numerically as far as we have looked.

This repository publishes the **code and the method**, plus a fixed-size **sample** of the output
under [`public/`](public/) — enough to reproduce most plots with no download at all. The full
tables are not published: they are gigabytes, and the code regenerates them in a few hours.

## The open questions

- **Conjecture.** `E(n,p) >= E(n,1/2)` for all `p` — that `p = 1/2` is the *global* minimum, not
  just a local one. [`RESEARCH_LOG.md`](RESEARCH_LOG.md) section 4 has the proof strategy that
  currently looks most promising, and the several that failed.
- **The sign of `F3`** at cusp points, and why the exceptions cluster where they do.
- Several structural facts that are verified numerically but not yet proved.

## Getting started

```bash
python3 -m venv .venv                       # use a native arm64 Python on Apple Silicon
.venv/bin/pip install numpy mpmath numba pyarrow matplotlib scipy

.venv/bin/python cusps_fast.py --nmax 1000 --workers 8 --out cusps/   # certified cusp tables
.venv/bin/python cusps_fast.py --merge --out cusps/

.venv/bin/python dump_ties.py --n 1000 --data data/ --workers 8       # every tie point, for plots
.venv/bin/python plotting/cusp_indicator.py --n 1000

.venv/bin/python star_check.py --nmax 5000 --workers 8 --out analysis/  # (★) at every switch point
```

`n = 1000` takes under two minutes on eight cores; `n <= 3000` takes about 2.5 hours.

## Where things are

| file | |
|---|---|
| [`RESEARCH_LOG.md`](RESEARCH_LOG.md) | **the substance** — definitions, proved facts, numerical results, proof approaches (including the failed ones), open questions, and a dated log of every run |
| [`README_cusps.md`](README_cusps.md) | how to run the pipeline, timings, output columns |
| [`binom_core.py`](binom_core.py) | the only implementation of the mathematics; everything else imports it |
| [`cusps_fast.py`](cusps_fast.py) | certified generator (numba screen + interval arithmetic) |
| [`dump_ties.py`](dump_ties.py) | Parquet export of every tie point |
| [`cusps_data.py`](cusps_data.py) | reader; derives the slope decomposition, gaps and indicators |
| [`CLAUDE.md`](CLAUDE.md) | working conventions, including the data publication policy |

Generated data is **not** committed — it is large and reproducible, and every ignored path has its
rebuild command in [`.gitignore`](.gitignore).

## License

MIT — see [`LICENSE`](LICENSE). The sample tables under `public/` are covered by the same terms.
