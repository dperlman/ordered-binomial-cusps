#!/bin/zsh
# Overnight run to n<=5000: generate (resumes from n<=3000), merge, gzip.
# Started 2026-09-21.  Log: run5000.log   Done marker: run5000.done
cd /Users/omgoleus/git/ordered-binomial-cusps
P=.venv/bin/python
echo "=== START $(date) ==="
echo "--- generate n=3..5000 (resumable; n<=3000 already present) ---"
$P cusps_fast.py --nmax 5000 --workers 8 --out cusps/ || { echo "GENERATE FAILED"; echo fail > run5000.done; exit 1; }
echo "--- merge $(date) ---"
$P cusps_fast.py --merge --out cusps/ || { echo "MERGE FAILED"; echo fail > run5000.done; exit 1; }
echo "--- gzip $(date) ---"
gzip -c cusps/cusps_all.csv > cusps_n5000.csv.gz || { echo "GZIP FAILED"; echo fail > run5000.done; exit 1; }
echo "--- summary ---"
ls -la cusps_n5000.csv.gz
wc -l cusps/cusps_all.csv
echo "interval checks total: $(wc -l < cusps/interval_checks.log)"
echo "UNRESOLVED rows: $(grep -c UNRESOLVED cusps/cusps_all.csv || true)"
echo "=== DONE $(date) ==="
echo ok > run5000.done
