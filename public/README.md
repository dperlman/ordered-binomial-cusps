# Public data tables

Built by `make_public_data.py`; do not edit by hand.

This repository publishes the **code and the method**, plus a **fixed-size sample** of the output.
It does not publish the full tables: they are gigabytes and they are reproducible from the code in
a few hours.  See `README_cusps.md` to generate them yourself.

| file | kind | what |
|---|---|---|
| `per_n_summary.csv` | result | one row per n: cusp counts, F3<0 counts, extremes, gap to E(1/2) |
| `cusps_decade.csv` | sample | every cusp for n = 100, 200, 500, 1000, 2000, 3000, 4000, 5000 |
| `cusps_F3_negative.csv` | sample | every F3<0 cusp with nearest-cusp metrics, n <= 3000 |

The samples are **pinned**: they stay at these n even as the local tables are extended, so the
repository does not grow with the frontier.  The result table is one row per n and may cover the
whole computed range.

Columns of the cusp tables:
`n,i,j,pstar,E,F3,F3_sign,S_minus,S_plus,slope_left,slope_right,certified_by`

For reference, the full local tables at the time of writing reach n=5000 with 4,421,154 cusps.
