import random

import pytest

from k3fq.labels import (poly_to_label, label_to_poly, parse_label, label_sort_key,
                         signed_cremona_letter_code, signed_class_to_int, cyclotomic_index,
                         balanced_code, balanced_value, transcendental_code, transcendental_coefficients,
                         cyclotomic_orders, num_cyclotomic_parts, cyclotomic_product,
                         split_cyclotomic, euler_phi, poly_mul, CYCLOTOMIC_ORDERS)


def test_codes():
    for m in list(range(-1000, 1000)) + [10**12, -10**12]:
        assert signed_class_to_int(signed_cremona_letter_code(m)) == m
        c = balanced_code(m)
        assert balanced_value(c) == m and (len(c) == 1 or c[0] != "n")
        assert balanced_value(balanced_code(m, 5)) == m
    # one letter for -13, ..., 12 and two for -351, ..., 324
    assert [len(balanced_code(m)) for m in (-13, 12, -14, 13, -351, 324, -352, 325)] == [1, 1, 2, 2, 2, 2, 3, 3]
    # with a common width, string order is numeric order
    vals = list(range(-351, 325))
    assert sorted(vals, key=lambda m: balanced_code(m, 2)) == vals


def test_transcendental_code():
    assert transcendental_code([0, 1, -1, 12, -13]) == "nomza"
    assert transcendental_code([1, 13]) == "nooa"
    assert transcendental_coefficients("nooa", 2) == [1, 13]
    for bad, k in [("nono", 2), ("nnnnnono", 2), ("nom", 2), ("", 1)]:
        with pytest.raises(ValueError):
            transcendental_coefficients(bad, k)


def test_cyclotomic_index():
    assert len(CYCLOTOMIC_ORDERS) == 43
    assert [num_cyclotomic_parts(k) for k in (0, 2, 4, 22)] == [1, 6, 24, 60280]
    for k in range(23):
        orders = [cyclotomic_orders(k, i) for i in range(num_cyclotomic_parts(k))]
        assert all(sum(euler_phi(m) for m in o) == k for o in orders)
        assert len({tuple(o) for o in orders}) == len(orders)
        assert [cyclotomic_index(o) for o in orders] == list(range(len(orders)))
        # exponent vectors (e_1, e_2, ...) decrease lexicographically
        vecs = [[o.count(m) for m in CYCLOTOMIC_ORDERS] for o in orders]
        assert vecs == sorted(vecs, reverse=True)
    assert cyclotomic_orders(22, 0) == [1] * 22
    assert cyclotomic_orders(22, 1) == [1] * 21 + [2]


def random_transcendental(rng, q, k):
    # not necessarily with roots on the unit circle; the labels don't care
    lower = [rng.randint(-10**6, 10**6) for _ in range(k)]
    return [q] + lower + lower[:-1][::-1] + [q] if k else [q]


def test_roundtrip():
    rng = random.Random(1)
    for _ in range(300):
        q = rng.choice([2, 3, 4, 5, 7, 8, 9, 25, 101])
        rho = rng.choice(range(0, 23, 2))
        orders = cyclotomic_orders(rho, rng.randrange(num_cyclotomic_parts(rho)))
        trans = random_transcendental(rng, q, (22 - rho) // 2)
        coeffs = poly_mul(cyclotomic_product(orders), trans)
        label = poly_to_label(coeffs)
        assert label.split(".")[:2] == [str(q), str(rho)]
        assert label_to_poly(label) == coeffs
        assert parse_label(label) == (q, orders, trans)
        # random transcendental parts this large have no cyclotomic factors
        assert split_cyclotomic(coeffs) == (orders, trans)


def test_examples():
    # 2(1 - T)^22 and 2(1 - T)^21(1 + T)
    assert poly_to_label(poly_mul(cyclotomic_product([1] * 22), [2])) == "2.22.a"
    assert poly_to_label(poly_mul(cyclotomic_product([1] * 21 + [2]), [2])) == "2.22.b"
    # 2 + 2T^22 = 2 Phi_4 Phi_44 is supersingular
    assert poly_to_label([2] + [0] * 21 + [2]) == "2.22.dilx"
    assert label_to_poly("2.22.dilx") == [2] + [0] * 21 + [2]
    # (1 - T)^16 Phi_4 (4 + T + 6T^2 + T^3 + 4T^4)
    P = poly_mul(cyclotomic_product([1] * 16 + [4]), [4, 1, 6, 1, 4])
    assert poly_to_label(P) == "4.18.e.ot"
    # (1 - T)^20 (2 + T + 2T^2)
    assert poly_to_label(poly_mul(cyclotomic_product([1] * 20), [2, 1, 2])) == "2.20.a.o"
    # b_2 = 17 needs two letters, so b_1 = -12 gets two as well
    P = poly_mul(cyclotomic_product([1] * 18), poly_mul([2, -3, 2], [2, -3, 2]))
    assert poly_to_label(P) == "4.18.a.nboe"
    assert label_to_poly("4.18.a.nboe") == P


@pytest.mark.parametrize("bad", [
    "2.p.a_a_a_a_a_a_a_a_a_a_a",   # old format
    "2.22",
    "02.22.a",
    "2.022.a",
    "1.22.a",
    "2.23.a",
    "2.24.a",
    "2.22.ab",                     # leading a in the cyclotomic code
    "2.22.dmcd",                   # index out of range
    "2.22.a.o",                    # no transcendental part when rho = 22
    "2.20.a",
    "2.20.a.",
    "2.20.a.o_o",                  # no separators
    "2.18.a.o_o",
    "2.18.a.ooo",                  # 3 letters for 2 coefficients
    "2.18.a.nono",                 # padded more than necessary
    "2.20.a.no",
    "2.20.A.o",
    "2.p.b_a_a_a_a_a_a_a_a_a_a",
])
def test_bad_labels(bad):
    with pytest.raises(ValueError):
        parse_label(bad)


def test_valid_edge_cases():
    assert parse_label("2.20.a.m")[2] == [2, -1, 2]
    assert parse_label("2.20.a.zz")[2] == [2, 324, 2]
    assert parse_label("2.18.a.naob")[2] == [2, -13, 14, -13, 2]
    assert parse_label("2.0.a.nnnnnnnnnnn") == (2, [], [2] + [0] * 21 + [2])


def test_sort_key():
    labels = ["2.22.b", "2.20.a.o", "2.20.ba.m", "2.20.c.m", "2.20.c.zz", "2.20.c.o", "3.0.a.nnnnnnnnnnn",
              "2.22.ba", "2.22.c", "2.20.c.ma"]
    assert sorted(labels, key=label_sort_key) == [
        "2.20.a.o", "2.20.c.ma", "2.20.c.m", "2.20.c.o", "2.20.c.zz", "2.20.ba.m", "2.22.b", "2.22.c", "2.22.ba",
        "3.0.a.nnnnnnnnnnn"]
