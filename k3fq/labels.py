"""
Labels for K3-type polynomials.

Write P = C * P_T, where C = (1 - T)^{e_1} prod_{m > 1} Phi_m^{e_m} is the
cyclotomic part, of degree rho (the geometric Picard rank), and P_T is the
transcendental part, which is palindromic of degree d = 22 - rho with
constant term q.  A label has the form

    q.rho.c.t      (k = d/2)

where

- ``c`` is the letter code of the index of the exponent vector
  (e_1, e_2, e_3, ...) among those of total degree rho, in decreasing
  lexicographic order (so ``a`` is (1 - T)^rho);
- ``t`` encodes the coefficients b_1, ..., b_k of T, ..., T^k in P_T: each is
  written in balanced base 26, with the letters a, ..., z standing for the
  digits -13, ..., 12 (so n = 0, o = 1, m = -1), padded with leading n's to
  a common width w, the smallest that fits all of them.  So ``t`` has
  length k * w, and when every -13 <= b_i <= 12 it has one letter per
  coefficient.

When rho = 22 there is no transcendental part and the label is ``q.22.c``.
The sign of the functional equation is (-1)^{e_1}, so it is not needed.

This module is pure Python, so that the LMFDB front end can share the logic.
"""
import re
from functools import lru_cache

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


def balanced_code(n, width=1):
    """
    ``n`` in balanced base 26, with the letters a-z standing for the digits
    -13, ..., 12, padded with leading n's (the digit 0) to ``width`` letters.

    >>> [balanced_code(m) for m in (0, 1, -1, 12, -13, 13, -14, 324, -351, 325)]
    ['n', 'o', 'm', 'z', 'a', 'oa', 'mz', 'zz', 'aa', 'oaa']
    >>> balanced_code(-1, 3)
    'nnm'
    """
    digits = []
    while n or not digits:
        r = (n + 13) % 26 - 13
        digits.append(r)
        n = (n - r) // 26
    digits.extend([0] * (width - len(digits)))
    return "".join(chr(110 + d) for d in reversed(digits))


def balanced_value(code):
    """
    Inverse of :func:`balanced_code`.

    >>> [balanced_value(c) for c in ('n', 'm', 'mz', 'nnm')]
    [0, -1, -14, -1]
    """
    n = 0
    for c in code:
        if not "a" <= c <= "z":
            raise ValueError("invalid letter code %r" % code)
        n = 26 * n + ord(c) - 110
    return n


def transcendental_code(lower):
    """
    The coefficients ``lower`` = [b_1, ..., b_k] in balanced base 26, all
    padded to the smallest common width.

    >>> transcendental_code([1, 6]), transcendental_code([-6, 12, -19, 26])
    ('ot', 'nhnzmuon')
    """
    width = max(len(balanced_code(b)) for b in lower)
    return "".join(balanced_code(b, width) for b in lower)


def transcendental_coefficients(code, k):
    """
    Inverse of :func:`transcendental_code`, for a code of ``k`` coefficients.
    Raises a ValueError if ``code`` is not in canonical form.
    """
    if not code or len(code) % k:
        raise ValueError("the transcendental code must have a multiple of %s letters" % k)
    width = len(code) // k
    chunks = [code[i:i + width] for i in range(0, len(code), width)]
    if width > 1 and all(c[0] == "n" for c in chunks):
        raise ValueError("the transcendental code is padded more than necessary")
    return [balanced_value(c) for c in chunks]


# ----------------------------------------------------------------------
# Integer polynomials, as coefficient lists with the constant term first

def _strip(f):
    f = list(f)
    while f and f[-1] == 0:
        f.pop()
    return f


def poly_mul(f, g):
    if not f or not g:
        return []
    h = [0] * (len(f) + len(g) - 1)
    for i, a in enumerate(f):
        if a:
            for j, b in enumerate(g):
                h[i + j] += a * b
    return h


def poly_divide(f, g):
    """
    The quotient f / g if g divides f exactly, else None.  The leading
    coefficient of g must be 1 or -1.
    """
    f, g = _strip(f), _strip(g)
    lead = g[-1]
    if lead not in (1, -1):
        raise ValueError("divisor must have leading coefficient 1 or -1")
    if len(f) < len(g):
        return [] if not f else None
    f = list(f)
    dg = len(g) - 1
    quo = [0] * (len(f) - dg)
    for i in range(len(f) - 1, dg - 1, -1):
        c = f[i] * lead
        if c:
            quo[i - dg] = c
            for j in range(dg + 1):
                f[i - dg + j] -= c * g[j]
    if any(f[:dg]):
        return None
    return quo


def _poly_eval(f, x):
    v = 0
    for c in reversed(f):
        v = v * x + c
    return v


def euler_phi(n):
    ans, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            ans -= ans // p
        p += 1
    if m > 1:
        ans -= ans // m
    return ans


# phi(m) >= sqrt(m/2), so m <= 2 * 22^2 suffices
CYCLOTOMIC_ORDERS = [m for m in range(1, 2 * DEGREE ** 2 + 1) if euler_phi(m) <= DEGREE]
_POSITION = {m: i for i, m in enumerate(CYCLOTOMIC_ORDERS)}


@lru_cache(maxsize=None)
def cyclotomic_factor(m):
    """
    Phi_m as a tuple of coefficients, normalized to have constant term 1
    (so 1 - T for m = 1).

    >>> cyclotomic_factor(1), cyclotomic_factor(6)
    ((1, -1), (1, -1, 1))
    """
    if m == 1:
        return (1, -1)
    f = [-1] + [0] * (m - 1) + [1]
    for d in range(1, m):
        if m % d == 0:
            f = poly_divide(f, cyclotomic_factor(d))
    if f[0] < 0:
        f = [-c for c in f]
    return tuple(f)


@lru_cache(maxsize=None)
def _tail_counts():
    """
    N[i][k] = the number of exponent vectors of degree k supported on
    CYCLOTOMIC_ORDERS[i:].
    """
    n = len(CYCLOTOMIC_ORDERS)
    N = [[0] * (DEGREE + 1) for _ in range(n + 1)]
    N[n][0] = 1
    for i in range(n - 1, -1, -1):
        ph = euler_phi(CYCLOTOMIC_ORDERS[i])
        for k in range(DEGREE + 1):
            N[i][k] = sum(N[i + 1][k - e * ph] for e in range(k // ph + 1))
    return N


def num_cyclotomic_parts(rho):
    """
    The number of cyclotomic parts of degree ``rho``.

    >>> [num_cyclotomic_parts(k) for k in (0, 2, 22)]
    [1, 6, 60280]
    """
    return _tail_counts()[0][rho]


def cyclotomic_index(orders):
    """
    The index of the cyclotomic part prod Phi_m (m in ``orders``, with
    multiplicity) among those of the same degree, ordered so that the
    exponent vectors (e_1, e_2, e_3, ...) decrease lexicographically.

    >>> [cyclotomic_index(o) for o in ([1, 1], [1, 2], [2, 2], [3], [4], [6])]
    [0, 1, 2, 3, 4, 5]
    """
    exps = {}
    for m in orders:
        if m not in _POSITION:
            raise ValueError("phi(%s) is larger than %s" % (m, DEGREE))
        exps[m] = exps.get(m, 0) + 1
    k = sum(euler_phi(m) for m in orders)
    if k > DEGREE:
        raise ValueError("cyclotomic part of degree %s > %s" % (k, DEGREE))
    N = _tail_counts()
    idx = 0
    for i, m in enumerate(CYCLOTOMIC_ORDERS):
        ph = euler_phi(m)
        e = exps.get(m, 0)
        # vectors agreeing so far but with a larger exponent of Phi_m come first
        for f in range(e + 1, k // ph + 1):
            idx += N[i + 1][k - f * ph]
        k -= e * ph
    return idx


def cyclotomic_orders(rho, idx):
    """
    Inverse of :func:`cyclotomic_index`: the sorted list of m (with
    multiplicity) for the cyclotomic part of degree ``rho`` with index ``idx``.

    >>> cyclotomic_orders(2, 4)
    [4]
    """
    N = _tail_counts()
    if not 0 <= rho <= DEGREE or not 0 <= idx < N[0][rho]:
        raise ValueError("no cyclotomic part of degree %s with index %s" % (rho, idx))
    orders = []
    k = rho
    for i, m in enumerate(CYCLOTOMIC_ORDERS):
        ph = euler_phi(m)
        for e in range(k // ph, -1, -1):
            n = N[i + 1][k - e * ph]
            if idx < n:
                break
            idx -= n
        orders.extend([m] * e)
        k -= e * ph
    return orders


def cyclotomic_product(orders):
    """
    prod Phi_m over ``orders``, normalized to have constant term 1.
    """
    C = [1]
    for m in orders:
        C = poly_mul(C, cyclotomic_factor(m))
    return C


def split_cyclotomic(coeffs):
    """
    Split the polynomial with coefficients ``coeffs`` as C * P_T, where C is
    a product of cyclotomic polynomials (normalized to have constant term 1)
    and P_T has no roots of unity among its roots.  Returns the sorted list
    of orders m (with multiplicity) and the coefficients of P_T.
    """
    f = _strip([int(c) for c in coeffs])
    v2, v3 = _poly_eval(f, 2), _poly_eval(f, 3)
    orders = []
    for m in CYCLOTOMIC_ORDERS:
        g = cyclotomic_factor(m)
        g2, g3 = _poly_eval(g, 2), _poly_eval(g, 3)
        while len(f) > len(g) - 1 and v2 % g2 == 0 and v3 % g3 == 0:
            h = poly_divide(f, g)
            if h is None:
                break
            f, v2, v3 = h, v2 // g2, v3 // g3
            orders.append(m)
    return orders, f


# ----------------------------------------------------------------------
# Labels

LABEL_RE = re.compile(r"^([1-9]\d*)\.(0|[1-9]\d*)\.([a-z]+)(?:\.([a-z]+))?$")


def poly_to_label(coeffs):
    """
    The label of the K3-type polynomial with coefficient list ``coeffs``
    (constant term first, length 23).

    >>> poly_to_label(poly_mul(cyclotomic_product([1] * 20), [2, 1, 2]))
    '2.20.a.o'
    >>> poly_to_label([2] + [0]*21 + [2])  # 2 Phi_4 Phi_44
    '2.22.dilx'
    """
    coeffs = [int(c) for c in coeffs]
    if len(coeffs) != DEGREE + 1:
        raise ValueError("expected %s coefficients" % (DEGREE + 1))
    q = coeffs[0]
    if q < 2 or coeffs[-1] not in (q, -q):
        raise ValueError("constant term must be q > 1 and leading coefficient +-q")
    orders, trans = split_cyclotomic(coeffs)
    return parts_to_label(q, orders, trans)


def parts_to_label(q, orders, trans):
    """
    The label of the class with cyclotomic orders ``orders`` (with
    multiplicity) and transcendental part with coefficients ``trans``
    (palindromic, constant term q), without multiplying them out.

    >>> parts_to_label(2, [1] * 20, [2, 1, 2])
    '2.20.a.o'
    """
    orders = sorted(orders)
    trans = [int(c) for c in trans]
    rho = sum(euler_phi(m) for m in orders)
    d = DEGREE - rho
    if len(trans) != d + 1 or trans != trans[::-1] or trans[0] != q:
        raise ValueError("transcendental part must be palindromic with constant term q")
    label = "%s.%s.%s" % (q, rho, cremona_letter_code(cyclotomic_index(orders)))
    if d:
        label += "." + transcendental_code(trans[1:d // 2 + 1])
    return label


def parse_label(label):
    """
    Returns ``(q, orders, trans)``: the cardinality of the base field, the
    sorted list of cyclotomic orders and the coefficients of the
    transcendental part.  Raises a ValueError for labels that are not in
    canonical form.
    """
    match = LABEL_RE.fullmatch(label)
    if not match:
        raise ValueError("labels have the form q.rho.c.t")
    q, rho, c, codes = match.groups()
    q, rho = int(q), int(rho)
    if q < 2:
        raise ValueError("q must be at least 2")
    if rho > DEGREE or rho % 2:
        raise ValueError("the geometric Picard rank must be even and at most %s" % DEGREE)
    idx = class_to_int(c)
    if cremona_letter_code(idx) != c:
        raise ValueError("the cyclotomic code %s has a leading a" % c)
    if idx >= num_cyclotomic_parts(rho):
        raise ValueError("there are only %s cyclotomic parts of degree %s" % (num_cyclotomic_parts(rho), rho))
    k = (DEGREE - rho) // 2
    if not k:
        if codes:
            raise ValueError("when rho = %s there is no transcendental part" % DEGREE)
        lower = []
    elif not codes:
        raise ValueError("geometric Picard rank %s needs a transcendental part" % rho)
    else:
        lower = transcendental_coefficients(codes, k)
    trans = [q] + lower + lower[:-1][::-1] + [q] if k else [q]
    return q, cyclotomic_orders(rho, idx), trans


def label_to_poly(label):
    """
    The coefficient list of the polynomial with the given label.

    >>> label_to_poly('2.22.b')[:3], label_to_poly('2.22.b')[-3:]
    ([2, -40, 378], [-378, 40, -2])
    """
    q, orders, trans = parse_label(label)
    return poly_mul(cyclotomic_product(orders), trans)


def label_sort_key(label):
    """
    A key for sorting labels: by q, geometric Picard rank, cyclotomic index
    and then the transcendental coefficients.
    """
    parts = label.split(".")
    k = (DEGREE - int(parts[1])) // 2
    lower = transcendental_coefficients(parts[3], k) if k else []
    return (int(parts[0]), int(parts[1]), class_to_int(parts[2])) + tuple(lower)
