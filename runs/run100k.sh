#!/bin/zsh
cd /Users/omgoleus/git/ordered-binomial-cusps
echo "=== START $(date) ==="
.venv/bin/python screen_collisions.py --nmax 100000 --workers 8 --save simplifying_n100000.csv
code=$?
echo "=== DONE $(date)  exit=$code ==="
ls -la simplifying_n100000.csv 2>/dev/null
echo $code > run100k.done
