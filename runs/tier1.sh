#!/bin/zsh
cd /Users/omgoleus/git/ordered-binomial-cusps
echo "=== [1/2] check_collisions --exhaustive 10000 (fixed: every pair within tol)  START $(date) ==="
.venv/bin/python check_collisions.py --exhaustive 10000 --workers 8
c1=$?
echo "=== [1/2] DONE $(date) exit=$c1 ==="
echo "=== [2/2] screen_collisions --no-fact-c --nmax 10000 (float-free Tier 1)  START $(date) ==="
.venv/bin/python screen_collisions.py --nmax 10000 --workers 8 --no-fact-c --save simplifying_n10000_tier1.csv
c2=$?
echo "=== [2/2] DONE $(date) exit=$c2 ==="
echo "$c1 $c2" > tier1.done
