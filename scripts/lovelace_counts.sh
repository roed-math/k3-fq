#!/bin/sh
# Count K3-type polynomials by geometric Picard rank, to choose the cutoffs
# (q, minimum geometric Picard rank) for the LMFDB tables.  Nothing but the
# counts is written: the enumeration output goes to /dev/null.
#
# Run on lovelace from the repository root:
#     git submodule update --init && make && sh scripts/lovelace_counts.sh
# The only expensive part is the full count over F_2 (about 5.5 million
# polynomials, all degrees up to 22); the rest takes about a minute.
set -e
J=${J:-32}
mkdir -p counts
sage -python k3fq.py counts -j "$J" 2 --max-degree 22 --output counts/q2_all.md
sage -python k3fq.py counts -j "$J" 3-64 --max-degree 6 --output counts/q3-64.md
