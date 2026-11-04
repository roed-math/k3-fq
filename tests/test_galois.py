import os

from sage.all import pari, ZZ

from k3fq.cyclotomic import R
from k3fq.galois import galois_label, chebyshev_real_poly, frobenius_kernel_weight
from k3fq.invariants import K3Poly

DATA = os.path.join(os.path.dirname(__file__), "data", "q2_sample.txt")


def test_small_degree_against_polgalois():
    # For degree <= 10 galois_label uses polgalois; check the wreath
    # certificate agrees by forcing it on degree 8 and 10 factors.
    from k3fq.galois import certify_full_wreath, real_galois_group, wreath_labels
    count = 0
    with open(DATA) as F:
        for line in F:
            kp = K3Poly([int(c) for c in line.split()])
            for f, e in kp.transcendental_factors():
                if f.degree() in (8, 10):
                    h = R(chebyshev_real_poly(f.list()))
                    H = real_galois_group(h)
                    label, full = galois_label(f)
                    if certify_full_wreath(f, h, H) is not None:
                        n, i = wreath_labels()[H]
                        assert label == "%sT%s" % (n, i)
                        assert full
                        count += 1
    assert count > 0


def test_kernel_weight():
    # pi = identity on 3 blocks, one block swapped: f has factors 2,1,1,1,1
    assert frobenius_kernel_weight([1, 1, 1], [2, 1, 1, 1, 1]) == 1
    # pi a 3-cycle that swaps: one factor of degree 6; sigma^3 = (1,1,1)
    assert frobenius_kernel_weight([3], [6]) == 3
    # pi = (2-cycle, swapping)(fixed, not swapping): m = 2, m/2 = 1 odd -> weight 2
    assert frobenius_kernel_weight([2, 1], [4, 1, 1]) == 2


def test_chebyshev():
    f = R([2, 3, 5, 3, 2])
    h = R(chebyshev_real_poly(f.list()))
    from sage.all import QQ, PolynomialRing
    t = PolynomialRing(QQ, "t").fraction_field().gen()
    assert f(t) == t ** 2 * h(t + 1 / t)
