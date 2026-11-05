"""
Tests for k3fq.invariants.  Run with ``sage -python -m pytest tests``.
"""
import pytest
from sage.all import ZZ, QQ, ComplexField, cyclotomic_polynomial, euler_phi, moebius

from k3fq.cyclotomic import R, T
from k3fq.invariants import K3Poly
from k3fq.labels import label_to_poly, poly_to_label
from k3fq.enumerate import times_one_minus_T2

CC = ComplexField(300)


def k3(P):
    return K3Poly(R(P).list()).compute()


def roots_with_mult(P):
    """
    Normalized Frobenius eigenvalues alpha_i (reciprocal roots), numerically.
    """
    ans = []
    for f, e in R(P).factor():
        if f.degree() > 0:
            ans.extend([1 / z for z in f.roots(CC, multiplicities=False)] * e)
    assert len(ans) == 22
    return ans


# A few test polynomials: transcendental parts over F_2 times cyclotomic parts
S4 = R([2, 2, 3, 2, 2])             # degree 4, from root-unitary k3f2-lines.txt
S2 = R([2, 3, 2])
EXAMPLES = [
    2 * (1 - T) ** 22,                                        # supersingular, rho = 22
    (1 - T) ** 2 * R(cyclotomic_polynomial(3)) ** 8 * S2 * (1 + T) ** 2,
    (1 - T) * (1 + T) * (1 - T) ** 16 * S4,                   # sign -1
    (1 - T) ** 2 * R(cyclotomic_polynomial(5)) ** 2 * R(cyclotomic_polynomial(12)) ** 2 * S2 * (1 + T) ** 2,
]


@pytest.fixture(scope="module", params=range(len(EXAMPLES)))
def example(request):
    P = EXAMPLES[request.param]
    assert P.degree() == 22 and P[0] == 2
    return P, k3(P)


def test_label_roundtrip(example):
    P, D = example
    assert label_to_poly(D["label"]) == P.list()
    assert poly_to_label(P.list()) == D["label"]


def test_sign(example):
    P, D = example
    assert D["sign"] == (-1) ** D["picard_rank"]
    assert P[22] == D["sign"] * 2


def test_point_counts_match_roots(example):
    P, D = example
    alphas = roots_with_mult(P)
    q = 2
    for n in range(1, 11):
        s = sum(a ** n for a in alphas)
        assert abs(s.imag()) < 1e-50
        count = 1 + q ** (2 * n) + q ** n * s.real()
        assert abs(count - D["k3_counts"][n - 1]) < 1e-40
        cubic = 1 + q ** n + q ** (2 * n) + q ** (3 * n) + q ** (4 * n) + q ** (2 * n) * s.real()
        assert abs(cubic - D["cubic_counts"][n - 1]) < 1e-30


def test_cyclotomic_part(example):
    P, D = example
    C = R(1)
    for m in D["cyclotomic_factors"]:
        C *= R(cyclotomic_polynomial(m))
    assert C.degree() == D["geometric_picard_rank"]
    S = R(D["transcendental_poly"])
    assert C * S in (P, -P)
    assert D["picard_rank"] == D["cyclotomic_factors"].count(1)


def test_picard_ranks(example):
    P, D = example
    alphas = roots_with_mult(P)
    for n, rho in zip(D["extension_degrees"], D["picard_ranks"]):
        assert rho == sum(1 for a in alphas if abs(a ** n - 1) < 1e-50)


def test_artin_tate(example):
    P, D = example
    alphas = roots_with_mult(P)
    q = 2
    for n, at in zip(D["extension_degrees"], D["artin_tate_values"]):
        val = q ** n
        for a in alphas:
            if abs(a ** n - 1) > 1e-50:
                val *= 1 - a ** n
        assert abs(val - at) < 1e-30 * abs(at)
    rho = D["picard_rank"]
    if rho:
        assert ZZ(D["ns_disc_class"]).is_squarefree()
        assert (D["artin_tate_value"] * D["ns_disc_class"] * (-1) ** (rho - 1)).is_square()


def test_ej_condition_matches_kedlaya(example):
    """
    Kedlaya-Sutherland's ej_test (k3-scripts/ej-test.sage), applied as in
    their scripts to P/(1-T).
    """
    P, D = example
    if D["picard_rank"] == 0:
        return
    pol = P // (1 - T)
    while pol % (1 - T) == 0:
        pol //= 1 - T
    assert ZZ(pol[0].sign() * pol(-1)).is_square() == D["ej_condition"]


def test_real_poly(example):
    P, D = example
    h = R(D["real_poly"])
    # P(T) = T^k h(T + 1/T) (times (1 - T^2) when the sign is -1)
    k = h.degree()
    from sage.all import PolynomialRing
    Rt = PolynomialRing(QQ, "t").fraction_field()
    t = Rt.gen()
    lhs = Rt(P(t))
    rhs = t ** k * h(t + 1 / t)
    if D["sign"] == -1:
        rhs *= (1 - t ** 2)
    assert lhs == rhs


def test_angles(example):
    P, D = example
    assert len(D["angles"]) == 11 + (D["cyclotomic_factors"].count(1) + D["cyclotomic_factors"].count(2)) // 2
    assert all(0 <= a <= 1 for a in D["angles"])


def test_supersingular():
    D = k3(2 * (1 - T) ** 22)
    assert D["is_supersingular"] and D["height"] is None
    assert D["newton_elevation"] == 11 and D["newton_coelevation"] == 0
    assert D["angle_rank"] == 0 and D["angle_corank"] == 0
    # AT_1 = 2 cannot be (square) * 2^(2 sigma): not realizable over F_2
    assert D["artin_tate_values"] == [2]
    assert not D["ej_condition"]
    assert D["realizable_k3"] == -1


def test_height_and_elevation():
    # Build polynomials of each height h from the Newton polygon: we just check
    # consistency of the definitions on random enumerated data.
    import os
    path = os.path.join(os.path.dirname(__file__), "data", "q2_sample.txt")
    with open(path) as F:
        for line in F:
            D = K3Poly([int(c) for c in line.split()]).compute()
            if not D["has_k3_newton_polygon"]:
                assert D["newton_elevation"] is None and D["newton_coelevation"] is None
            elif not D["is_supersingular"]:
                h = D["height"]
                assert D["newton_elevation"] == h - 1
                assert D["newton_coelevation"] == 12 - h
            assert D["is_ordinary"] == (D["poly"][1] % 2 != 0)


def test_count_growth_definition():
    path = __import__("os").path.join(__import__("os").path.dirname(__file__), "data", "q2_sample.txt")
    with open(path) as F:
        for line in F:
            kp = K3Poly([int(c) for c in line.split()])
            D = kp.compute()
            N = kp.counts("k3", 20)
            closed = [sum(moebius(n // d) * N[d] for d in ZZ(n).divisors()) for n in range(1, 21)]
            assert D["has_count_growth"] == all(c >= 0 for c in closed)
            assert D["has_nonnegative_counts"] == all(c >= 0 for c in N[1:])


def test_sign_minus_construction():
    c = [2, 1, 0, 1, 2]
    assert R(times_one_minus_T2(c)) == (1 - T ** 2) * R(c)


def test_fourfold_counts(example):
    # the F_q counts of the fourfolds are affine in the K3 count
    P, D = example
    q, n = 2, D["k3_count"]
    assert D["cubic_count"] == 1 + q**2 + q**4 + q * n == D["cubic_counts"][0]
    assert D["gm_count"] == D["cubic_count"] + q**2 == D["gm_counts"][0]
    assert D["verra_count"] == 1 + q + 2 * q**2 + q**3 + q**4 + q * n == D["verra_counts"][0]
