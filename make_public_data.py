"""make_public_data.py -- build the small CSV tables that ARE committed to the repo.

    python make_public_data.py [--cusps cusps/cusps_all.csv] [--analysis analysis/]
                               [--out public/] [--dry-run] [--allow-large]

Writes:
    public/per_n_summary.csv     RESULT   one row per n, full local range
    public/cusps_decade.csv      SAMPLE   every cusp for n in DECADE (pinned)
    public/cusps_F3_negative.csv SAMPLE   every F3<0 cusp with nearest-cusp metrics, n <= F3_NMAX
    public/README.md             what these are

WHAT GOES IN GIT, AND WHY IT DOES NOT GROW.
The repo offers two things: the code and the accumulated knowledge to generate these tables at any
n, and a fixed-size sample of what they look like.  It deliberately does NOT offer the full tables:
those are gigabytes, they are reproducible from the code, and nobody ever downloaded them when they
were published as release assets.  So committed data is divided by SHAPE, not by size:

  RESULT  -- one row per n.  Aggregates and extremes; no raw rows.  Grows LINEARLY in n, ~164
             bytes per n, so it may track the local frontier freely: it IS the finding.
  SAMPLE  -- one row per cusp.  Raw rows, so it shows what the output actually looks like.  Grows
             as n^2 (or per n added to DECADE), so it is PINNED and does not move when the local
             frontier moves.  The pins below are final; raising one is a deliberate decision that
             must be re-measured against the size bands, not a routine consequence of a longer run.

SIZE BANDS (enforced below, not merely documented -- prose is advice a future agent can rationalise
past, an assert is not):
    <= 1 MB   fine
    1-2 MB    fine, mention it in the commit message
    2-5 MB    HUMAN DISCRETION -- refused unless --allow-large is passed by a person
    >  5 MB   HARD STOP -- not committed, at any size, for any reason

Always regenerate with this script; never hand-edit.  Line endings are forced to LF everywhere
(Python's csv module defaults to the 'excel' dialect, i.e. CRLF, whatever the platform).
"""
import argparse, csv, os, sys

DECADE  = (100, 200, 500, 1000, 2000, 3000, 4000, 5000)   # PINNED 2026-09-22 -- do not extend
F3_NMAX = 3000                                            # PINNED 2026-09-22 -- n^2 growth

SOFT_MB, DISCRETION_MB, HARD_MB = 1.0, 2.0, 5.0

def check_size(path, allow_large):
    mb = os.path.getsize(path)/1e6
    if mb > HARD_MB:
        sys.exit(f"REFUSED: {path} is {mb:.2f} MB, over the {HARD_MB} MB hard limit. "
                 f"It does not go in git at any size. Pin it lower or reshape it as a RESULT.")
    if mb > DISCRETION_MB and not allow_large:
        sys.exit(f"REFUSED: {path} is {mb:.2f} MB, in the {DISCRETION_MB}-{HARD_MB} MB discretionary "
                 f"band. A person must pass --allow-large; an agent must not decide this alone.")
    band = ("ok" if mb <= SOFT_MB else
            "note in commit" if mb <= DISCRETION_MB else "DISCRETIONARY (--allow-large)")
    return mb, band

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cusps", default="cusps/cusps_all.csv")
    ap.add_argument("--analysis", default="analysis")
    ap.add_argument("--out", default="public")
    ap.add_argument("--dry-run", action="store_true", help="report sizes, write nothing")
    ap.add_argument("--allow-large", action="store_true",
                    help="permit a file in the 2-5 MB discretionary band (a PERSON passes this)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    tmp = a.out if not a.dry_run else os.path.join(a.out, ".dryrun")
    os.makedirs(tmp, exist_ok=True)
    written = []

    # --- SAMPLE: every cusp for the pinned DECADE values of n -------------------
    dec = os.path.join(tmp, "cusps_decade.csv")
    with open(a.cusps) as fi, open(dec, "w", newline="") as fo:
        rd = csv.reader(fi); w = csv.writer(fo, lineterminator="\n")
        w.writerow(next(rd)); n_rows = 0; n_all = 0; n_max = 0
        for r in rd:
            n_all += 1; n = int(r[0])
            if n > n_max: n_max = n
            if n in DECADE: w.writerow(r); n_rows += 1
    written.append((dec, f"{n_rows} rows, n={','.join(map(str, DECADE))}", "SAMPLE (pinned)"))

    # --- RESULT: one row per n, whole local range -------------------------------
    src = os.path.join(a.analysis, "per_n_summary.csv")
    dst = os.path.join(tmp, "per_n_summary.csv")
    if os.path.exists(src):
        # copy, normalising line endings: files written before 2026-09-22 are CRLF (Python's csv
        # module defaults to the 'excel' dialect), and everything committed here must be LF.
        with open(src, newline='') as fi, open(dst, 'w', newline='') as fo:
            k = -1
            for line in fi: fo.write(line.replace('\r\n', '\n')); k += 1
        written.append((dst, f"{k} rows, n=3..{n_max}", "RESULT (tracks the frontier)"))
    else:
        sys.exit(f"MISSING {src} -- run analyze_cusps.py first")

    # --- SAMPLE: F3<0 cusps, pinned at F3_NMAX ----------------------------------
    src = os.path.join(a.analysis, "cusps_F3_negative.csv")
    dst = os.path.join(tmp, "cusps_F3_negative.csv")
    if os.path.exists(src):
        with open(src) as fi, open(dst, "w", newline="") as fo:
            rd = csv.reader(fi); w = csv.writer(fo, lineterminator="\n")
            w.writerow(next(rd)); k = 0
            for r in rd:
                if int(r[0]) <= F3_NMAX: w.writerow(r); k += 1
        written.append((dst, f"{k} rows, n<={F3_NMAX}", "SAMPLE (pinned)"))
    else:
        sys.exit(f"MISSING {src} -- run analyze_cusps.py first")

    print(f"{'file':28s} {'size':>9}  {'band':<22} what")
    for path, what, kind in written:
        mb, band = check_size(path, a.allow_large)
        print(f"{os.path.basename(path):28s} {mb*1000:7.0f} KB  {band:<22} {kind}: {what}")

    with open(os.path.join(tmp, "README.md"), "w") as fh:
        fh.write(f"""# Public data tables

Built by `make_public_data.py`; do not edit by hand.

This repository publishes the **code and the method**, plus a **fixed-size sample** of the output.
It does not publish the full tables: they are gigabytes and they are reproducible from the code in
a few hours.  See `README_cusps.md` to generate them yourself.

| file | kind | what |
|---|---|---|
| `per_n_summary.csv` | result | one row per n: cusp counts, F3<0 counts, extremes, gap to E(1/2) |
| `cusps_decade.csv` | sample | every cusp for n = {', '.join(map(str, DECADE))} |
| `cusps_F3_negative.csv` | sample | every F3<0 cusp with nearest-cusp metrics, n <= {F3_NMAX} |

The samples are **pinned**: they stay at these n even as the local tables are extended, so the
repository does not grow with the frontier.  The result table is one row per n and may cover the
whole computed range.

Columns of the cusp tables:
`n,i,j,pstar,E,F3,F3_sign,S_minus,S_plus,slope_left,slope_right,certified_by`

For reference, the full local tables at the time of writing reach n={n_max} with {n_all:,} cusps.
""")
    if a.dry_run:
        print(f"\n--dry-run: wrote nothing to {a.out}/ (scratch in {tmp}/)")
    else:
        print(f"-> {a.out}/")

if __name__ == "__main__":
    main()
