"""
Labels for K3-type polynomials.

A label has the form ``q.s.c1_c2_..._c11`` where ``s`` is ``p`` (sign +1) or
``m`` (sign -1) and ``ci`` is the signed Cremona letter code of the
coefficient a_i of P(T) = q + a_1 T + ... + a_22 T^22.  The upper half of the
coefficients is determined by a_{22-i} = sign * a_i.
"""

DEGREE = 22
HALF = DEGREE // 2


def cremona_letter_code(n):
    """
    Base 26 encoding of a nonnegative integer with digits a-z (a = 0).

    >>> [cremona_letter_code(n) for n in (0, 1, 25, 26, 27, 44)]
    ['a', 'b', 'z', 'ba', 'bb', 'bs']
    """
    if n < 0:
        raise ValueError("n must be nonnegative")
    if n == 0:
        return "a"
    digits = []
    while n:
        n, r = divmod(n, 26)
        digits.append(chr(97 + r))
    return "".join(reversed(digits))


def class_to_int(code):
    """
    Inverse of :func:`cremona_letter_code`.

    >>> class_to_int('bs')
    44
    """
    n = 0
    for c in code:
        if not "a" <= c <= "z":
            raise ValueError("invalid letter code %r" % code)
        n = 26 * n + ord(c) - 97
    return n


def signed_cremona_letter_code(m):
    """
    Signed variant used in LMFDB abelian variety labels: negative numbers
    get a leading ``a``.

    >>> [signed_cremona_letter_code(m) for m in (0, 3, -3, -44)]
    ['a', 'd', 'ad', 'abs']
    """
    if m >= 0:
        return cremona_letter_code(m)
    return "a" + cremona_letter_code(-m)


def signed_class_to_int(code):
    """
    >>> [signed_class_to_int(c) for c in ('a', 'd', 'ad', 'abs')]
    [0, 3, -3, -44]
    """
    if code == "a":
        return 0
    if code.startswith("a"):
        return -class_to_int(code[1:])
    return class_to_int(code)


def poly_to_label(coeffs):
    """
    The label of the K3-type polynomial with coefficient list ``coeffs``
    (constant term first, length 23).

    >>> poly_to_label([2] + [0]*21 + [2])
    '2.p.a_a_a_a_a_a_a_a_a_a_a'
    """
    if len(coeffs) != DEGREE + 1:
        raise ValueError("expected %s coefficients" % (DEGREE + 1))
    q = int(coeffs[0])
    lead = int(coeffs[-1])
    if lead == q:
        s = "p"
    elif lead == -q:
        s = "m"
    else:
        raise ValueError("leading coefficient must be +-q")
    codes = "_".join(signed_cremona_letter_code(int(c)) for c in coeffs[1:HALF + 1])
    return "%s.%s.%s" % (q, s, codes)


def label_to_poly(label):
    """
    The coefficient list of the polynomial with the given label.

    >>> label_to_poly('2.m.ab_a_a_a_a_a_a_a_a_a_a')[:3], label_to_poly('2.m.ab_a_a_a_a_a_a_a_a_a_a')[-3:]
    ([2, -1, 0], [0, 1, -2])
    """
    q, s, codes = label.split(".")
    q = int(q)
    sign = {"p": 1, "m": -1}[s]
    lower = [signed_class_to_int(c) for c in codes.split("_")]
    if len(lower) != HALF:
        raise ValueError("expected %s coefficient codes" % HALF)
    if sign == -1 and lower[-1] != 0:
        raise ValueError("middle coefficient must vanish when the sign is -1")
    coeffs = [q] + lower + [sign * c for c in reversed(lower[:-1])] + [sign * q]
    return coeffs
