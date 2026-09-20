"""
make_public_data.py -- build the small CSV tables that ARE committed to the repo (Tier 1/2).

    python make_public_data.py [--cusps cusps/cusps_all.csv] [--analysis analysis/] [--out public/]

Writes:
    public/per_n_summary.csv     one row per n (copied from analyze_cusps.py output)
    public/cusps_decade.csv      every cusp for n in DECADE
    public/negF3_neighbors.csv   every F3<0 cusp with its nearest-cusp metrics  (Tier 2)
    public/README.md             what these are and how to get the full tables

These are the only data files in git.  The full tables go to GitHub Releases -- see the data
publication policy in CLAUDE.md.  Always regenerate with this script; never hand-edit.
"""
import argparse, csv, os, shutil

DECADE = (100, 200, 500, 1000, 2000)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cusps", default="cusps/cusps_all.csv")
    ap.add_argument("--analysis", default="analysis")
    ap.add_argument("--out", default="public")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    with open(a.cusps) as fi, open(os.path.join(a.out, "cusps_decade.csv"), "w", newline="") as fo:
        rd = csv.reader(fi); w = csv.writer(fo)
        w.writerow(next(rd)); n_rows = 0
        for r in rd:
            if int(r[0]) in DECADE: w.writerow(r); n_rows += 1
    print(f"cusps_decade.csv: {n_rows} rows (n={','.join(map(str, DECADE))})")

    for name, tier in (("per_n_summary.csv", 1), ("negF3_neighbors.csv", 2)):
        src = os.path.join(a.analysis, name)
        if os.path.exists(src):
            shutil.copyfile(src, os.path.join(a.out, name))
            with open(src) as fh: k = sum(1 for _ in fh) - 1
            print(f"{name}: {k} rows (tier {tier})")
        else:
            print(f"{name}: MISSING ({src}) -- run analyze_cusps.py first")

    with open(os.path.join(a.out, "README.md"), "w") as fh:
        fh.write(f"""# Public data tables

Built by `make_public_data.py`; do not edit by hand.

| file | what |
|---|---|
| `per_n_summary.csv` | one row per n: cusp counts, F3<0 counts, extremes, gap to E(1/2) |
| `cusps_decade.csv` | every cusp for n = {', '.join(map(str, DECADE))} |
| `negF3_neighbors.csv` | every F3<0 cusp with nearest-cusp distance metrics |

Columns of the cusp tables:
`n,i,j,pstar,E,F3,F3_sign,S_minus,S_plus,slope_left,slope_right,certified_by`

These are extracts.  The full tables (n=3..2000, 707,417 cusps) are attached to the GitHub
Releases of this repository, not stored in git.  Everything else is regenerated locally --
see `README_cusps.md`.
""")
    print(f"-> {a.out}/")

if __name__ == "__main__":
    main()
