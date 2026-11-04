import random

from k3fq.labels import poly_to_label, label_to_poly, signed_cremona_letter_code, signed_class_to_int


def test_codes():
    for m in list(range(-1000, 1000)) + [10**12, -10**12]:
        assert signed_class_to_int(signed_cremona_letter_code(m)) == m


def test_roundtrip():
    rng = random.Random(1)
    for _ in range(200):
        q = rng.choice([2, 3, 4, 5, 7, 8, 9, 25, 101])
        sign = rng.choice([1, -1])
        lower = [rng.randint(-10**6, 10**6) for _ in range(11)]
        if sign == -1:
            lower[-1] = 0
        coeffs = [q] + lower + [sign * c for c in reversed(lower[:-1])] + [sign * q]
        assert label_to_poly(poly_to_label(coeffs)) == coeffs
