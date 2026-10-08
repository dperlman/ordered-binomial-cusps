#!/bin/zsh
cd /Users/omgoleus/git/ordered-binomial-cusps
echo "=== START $(date) ==="
.venv/bin/python verify_fact_c.py --nmin 39950 --nmax 40050 \
    --ns 27719 45359 50399 55439 65519 83159 98279 \
    --workers 8 --catalogue simplifying_n100000.csv
code=$?
echo "=== DONE $(date)  exit=$code ==="
echo $code > spotcheck.done
