# runs/

The launch scripts of the long batch runs whose results are cited in [FACTS.md](../FACTS.md) and
[RESEARCH_LOG.md](../RESEARCH_LOG.md): the exact invocations, kept as the record of how those
results were produced. Each `cd`s to the repository root and writes a `.log` (via the shell
redirect it was started with) and a `.done` marker holding its exit code; both are git-ignored.

| script | what it ran | results |
|---|---|---|
| `run5000.sh` | `cusps_fast.py` to n ≤ 5000 (resuming from n ≤ 3000), then merge and gzip (2026-09-21/22) | the complete certified cusp tables (FACTS S4, and the n ≤ 5000 rows S5–S10) |
| `tier1.sh` | `check_collisions.py --exhaustive 10000`, then the float-free `screen_collisions.py --no-fact-c --nmax 10000` | no tie-point collisions for n ≤ 10,000 by two independent methods (FACTS S1, S3) |
| `run100k.sh` | `screen_collisions.py --nmax 100000 --save simplifying_n100000.csv` (2026-09-24) | no collisions for n ≤ 100,000 under Fact C, and the reducing-point catalogue (FACTS S2, S3) |
| `spotcheck.sh` | `verify_fact_c.py` brute-force spot checks at n = 39,950–40,050 and seven larger n | the independent check of the n ≤ 100,000 screen (RESEARCH_LOG 2026-09-24) |

They run from the repository root (`zsh runs/<script>.sh`); the paths inside are relative to it.
