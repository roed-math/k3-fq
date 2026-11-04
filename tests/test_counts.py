"""
The counting identities in k3fq/counts.py, checked against Kedlaya's data
and against direct enumeration.
"""
from k3fq.counts import cyclotomic_counts, weil_count, transcendental_counts, picard_table
from k3fq.enumerate import cyclotomic_products

# T(2, 2k) and T(3, 2k) for k <= 5, from root-unitary/k3-scripts/k3f{2,3}-lines.txt
KEDLAYA_T2 = [1, 4, 18, 62, 200, 546]
KEDLAYA_T3 = [1, 8, 56, 318, 1620, 7230]


def test_cyclotomic_counts():
    C, Cplus = cyclotomic_counts()
    for e in range(0, 13, 2):
        prods = cyclotomic_products(e)
        assert len(prods) == C[e]
        assert sum(1 for f in prods if _phi1_mult(f) % 2 == 0) == Cplus[e]
    assert C[22] == 60280


def _phi1_mult(f):
    from k3fq.cyclotomic import T
    e = 0
    while f % (1 - T) == 0:
        f //= 1 - T
        e += 1
    return e


def test_transcendental_counts():
    C, Cplus = cyclotomic_counts()
    for q, expected in [(2, KEDLAYA_T2), (3, KEDLAYA_T3)]:
        W = [weil_count(2 * k, q, threads=2)[0] for k in range(len(expected))]
        assert transcendental_counts(W, Cplus) == expected
    for q in (4, 5, 7, 8, 9):
        W = [weil_count(2 * k, q, threads=2)[0] for k in range(2)]
        assert transcendental_counts(W, Cplus)[1] == 4 * q - 4


def test_picard_table_matches_enumeration():
    # enumerate_elevated(2, 14) wrote 1,546,496 polynomials, enumerate_elevated(3, 16) 3,413,708
    C, _ = cyclotomic_counts()
    assert picard_table([1, 4, 18, 62, 200], C)[14] == 1546496
    assert picard_table(KEDLAYA_T3[:4], C)[16] == 3413708
