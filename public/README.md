# Public data tables

Built by `make_public_data.py`; do not edit by hand.

| file | what |
|---|---|
| `per_n_summary.csv` | one row per n: cusp counts, F3<0 counts, extremes, gap to E(1/2) |
| `cusps_decade.csv` | every cusp for n = 100, 200, 500, 1000, 2000 |
| `negF3_neighbors.csv` | every F3<0 cusp with nearest-cusp distance metrics |

Columns of the cusp tables:
`n,i,j,pstar,E,F3,F3_sign,S_minus,S_plus,slope_left,slope_right,certified_by`

These are extracts.  The full tables (n=3..2000, 707,417 cusps) are attached to the GitHub
Releases of this repository, not stored in git.  Everything else is regenerated locally --
see `README_cusps.md`.
