# k3-fq: motives of K3 type over finite fields

Scripts for building the LMFDB tables of motives of K3 type over finite fields, up to isogeny.
Each row of the main table corresponds to a polynomial

    P(T) = q + a_1 T + ... + a_22 T^22  in Z[T],   a_22 = ±q,

with all roots on the unit circle.  If X is a K3 surface over F_q, then
P(T) = q det(1 - q^{-1} Frob T | H^2).  The table schemas are in [`schemas/`](schemas/):

* [`k3_fq_isog`](schemas/k3_fq_isog.md): one row per polynomial.
* [`k3_fq_lattices`](schemas/k3_fq_lattices.md): candidate geometric Néron–Severi lattices with their Frobenius isometry.  Not computed here yet.

The enumeration uses Kiran Kedlaya's [root-unitary](https://github.com/kedlaya/root-unitary) code directly (a git submodule in `external/`), not through Sage's `weil_polynomials`.  Its OpenMP driver `c-script/weilpoly.c` parallelizes the search tree itself.  The structure of the rest is loosely based on [abvar-fq](https://github.com/roed314/abvar-fq).  It is much simpler, because nothing is decomposed into simple factors: each polynomial is processed independently, apart from the twist and base change stages.

## Setup

    git clone --recursive <this repository>
    make            # builds bin/weilpoly against the FLINT shipped with Sage
    make test       # sage -python -m pytest tests

Everything else runs under `sage -python` (Sage 10.x).  The `fields` stage uses Magma for Galois groups that the cheap certificate cannot decide.  The `kummer` stage reads abelian surfaces from the public LMFDB mirror.

## Running

All stages take `--data DIR` (default `data/`) and `-j JOBS`.

    sage -python k3fq.py -j 14 enumerate 2                           # all of them: data/polys/q2_p.txt, q2_m.txt
    sage -python k3fq.py -j 14 enumerate 7 --min-geometric-rank 10   # only geometric Picard rank >= 10
    sage -python k3fq.py -j 14 basic data/polys/q2_*.txt            # invariants, in chunks (resumable)
    sage -python k3fq.py -j 14 fields 2                             # Galois groups (Magma fallback), polredbest
    sage -python k3fq.py -j 14 twists 2                             # twist classes
    sage -python k3fq.py basechange                                 # primitive models (all q present)
    sage -python k3fq.py kummer 3                                   # Kummer candidates (needs network)
    sage -python k3fq.py export 2                                   # data/export/k3_fq_isog.q2.txt

The export file is in psycodict's `copy_from`/`reload` format: a header line with column names, a line of types, a blank line, then `|`-separated rows.

### Stages

| stage | what it does | output |
| --- | --- | --- |
| `enumerate` | Sign +1: `weilpoly 22 1 q`.  Sign -1: `(1 - T^2)` times the output of `weilpoly 20 1 q`.  With `--min-geometric-rank R`, it multiplies the transcendental parts of degree `2k <= 22 - R` (output of `weilpoly 2k 1 q` with no cyclotomic factor) by all cyclotomic products of degree `22 - 2k`, as in Kedlaya–Sutherland | `polys/q{q}_*.txt` |
| `basic` | everything that depends only on `P`: label, real polynomial, angles, angle rank, Newton polygon (slopes, height, elevation), cyclotomic factors, Picard ranks over extensions, Artin–Tate values, the Elsenhans–Jahnel condition, transcendental factorization and (usually) its Galois group, point counts for the four variety types, count positivity/growth, and the realizability flags that follow from these necessary conditions | `basic/q{q}/*.txt` |
| `fields` | Galois groups the certificate missed (Magma's `GaloisGroup`); `polredbest` for each transcendental factor | `fields/q{q}.txt` |
| `twists` | geometric class (twist class), its size, and for classes of at most 50 members the minimal twist degrees with the common base change | `twists/q{q}.txt` |
| `basechange` | primitive models: base changes of the classes over every subfield present in the data directory | `basechange/q{q}.txt` |
| `kummer` | abelian surfaces over F_q (from `av_fq_isog`) whose wedge square has the same transcendental part | `kummer/q{q}.txt` |
| `export` | joins everything; columns from stages that have not been run are null | `export/k3_fq_isog.q{q}.txt` |

### Where to run

Generate the real data on lovelace.  mimir is shared and often heavily loaded, so use it only for development and small tests.

## Uploading to the LMFDB

[`test_upload/k3_fq_isog_test.txt`](test_upload/k3_fq_isog_test.txt) is a 12-row test file in the format read by psycodict's `copy_from`.  It contains the two supersingular classes over F_2 that the LMFDB tests (`lmfdb/k3fq/test_k3fq.py`) look up, plus classes with assorted features (both signs, Kummer candidates, a base change from F_2 to F_4 together with its model, a non-K3 Newton polygon, and a reducible transcendental part).  `scripts/make_test_upload.py` regenerates it from any export file.

To create the table and load a file, from an LMFDB checkout with write access (this has not been run against a real database):

```python
import sys
from collections import defaultdict
sys.path.append("/path/to/k3-fq")
from k3fq.columns import ISOG_COLUMNS
from lmfdb import db

search_columns = defaultdict(list)
for col, typ in ISOG_COLUMNS:
    search_columns[typ].append(col)
db.create_table("k3_fq_isog", dict(search_columns), label_col="label", sort=["q", "poly"],
                table_description="Motives of K3 type over finite fields, up to isogeny",
                force_description=False)
db.k3_fq_isog.copy_from("test_upload/k3_fq_isog_test.txt")

from lmfdb.k3fq.stats import K3FqStats
K3FqStats().setup()   # statistics for /Motive/K3/Fq/stats
```

For the full table you will also want indexes, e.g. on `(q, poly)` (the default sort) and on `poly_str` with `text_pattern_ops` (for initial-coefficient searches).

## Notes on the algorithms

* **Galois groups.**  A transcendental factor `f` of degree `2k` is palindromic with no real roots, so `f(T) = T^k h(T + 1/T)` and `Gal(f)` lies in `C2 wr Gal(h)`.  `Gal(h)` has degree at most 11, so PARI's `polgalois` handles it.  We certify that `Gal(f)` is the full wreath product by finding a prime whose Frobenius gives a suitable element of the kernel; see `k3fq/galois.py`.  The certificate reads the cycle structure off the factorizations of `f` and `h` modulo the prime.  The label of `C2 wr kTj` comes from `k3fq/data/wreath_labels.txt`, which `magma/wreath_labels.m` generates.  When the certificate holds, the angle rank is exactly `k`.  Otherwise we use PARI's `lindep` at high precision, as abvar-fq does.  In a sample of 60 factors of degree at least 12 over F_2, the certificate succeeded for 48, and Magma confirmed all 48 labels; the 12 failures had genuinely smaller Galois groups.
* **Point count checks.**  `has_nonnegative_counts` and `has_count_growth` are checked for all `n` with `q^n < 10^4`.  Beyond that, both hold automatically for every variety type: the top power of `q` dominates the trace term, which has absolute value at most `22 q^{wn}`.
* **Twists.**  Twist classes are found by hashing the multiset of `{M θ}` over the transcendental angles, with `M = lcm{m : φ(m) ≤ 44}`; see `k3fq/twists.py`.  The degrees over which the cyclotomic parts agree, and those over which the transcendental parts agree, are both closed under multiples.  So the minimal twist degrees are minimal lcms of the minimal degrees for the two parts.  These are found exactly for the cyclotomic part and numerically for the transcendental part.  Each degree with `q^n < 2^128` is then checked by exact base change, and its base-change label recorded.  Base change uses Newton's identities on the power sums `Tr(Frob^{kn})`.
* **Sign -1.**  This forces `(1 - T)(1 + T) | P`, and in fact the Picard rank over F_q is odd exactly when the sign is -1.

## Performance

On mimir (shared and heavily loaded at the time), the `basic` stage takes about 12–15 ms per polynomial on one core.  Most of that is the Galois certificate and root finding.

According to the discussion notes, over F_2 there are just under 4 million polynomials with sign +1 and about 1.25 million with sign -1 (I have not rerun the full enumeration here).  Most of them have large cyclotomic parts: there are 60,280 products of cyclotomic polynomials of degree 22, and 1,546,496 of the polynomials over F_2 have geometric Picard rank at least 14.

## Status

Not yet implemented, so exported as null:

* LMFDB number field labels (`number_fields`).  This needs `polredabs`, which is slow for large degree.
* `polarization_degrees` (Taelman's theorem for ordinary classes).
* The lattice table, and the lattice summary columns.
* Realizability beyond the necessary conditions.

Open questions are listed at the end of [`schemas/k3_fq_isog.md`](schemas/k3_fq_isog.md).
