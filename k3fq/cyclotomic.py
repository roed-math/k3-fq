"""
Bookkeeping for cyclotomic factors of K3-type polynomials.
"""
from functools import lru_cache
from math import gcd

from sage.all import ZZ, PolynomialRing, cyclotomic_polynomial, euler_phi

R = PolynomialRing(ZZ, "T")
T = R.gen()


def lcm(a, b):
    return a * b // gcd(a, b)


@lru_cache(maxsize=None)
def cyclotomic_orders_up_to(n):
    """
    The list of m with phi(m) <= n.
    """
    # phi(m) >= sqrt(m/2), so m <= 2 n^2 suffices
    return [m for m in range(1, 2 * n * n + 3) if euler_phi(m) <= n]


@lru_cache(maxsize=None)
def _order_table(n):
    """
    Map from coefficient tuples of Phi_m (normalized to positive leading
    coefficient) to m, for phi(m) <= n.
    """
    return {tuple(cyclotomic_polynomial(m).list()): m for m in cyclotomic_orders_up_to(n)}


def cyclotomic_order(f):
    """
    If ``f`` is (up to sign) a cyclotomic polynomial Phi_m, return m; else 0.
    """
    if f.leading_coefficient() < 0:
        f = -f
    return _order_table(max(f.degree(), 1)).get(tuple(f.list()), 0)


def phi_at_one(m):
    """
    Phi_m(1): p if m is a power of the prime p, 0 if m = 1, else 1.
    """
    if m == 1:
        return 0
    F = ZZ(m).factor()
    if len(F) == 1:
        return F[0][0]
    return 1


def split_cyclotomic(P):
    """
    Split ``P`` (with P(0) > 0) as ``C * S`` with ``C`` a product of
    cyclotomic polynomials normalized to have constant term 1 and ``S``
    without roots of unity.

    Returns ``(orders, C, S)`` where ``orders`` is the sorted list of m (with
    multiplicity) such that C = prod Phi_m up to sign.
    """
    C = P.cyclotomic_part()
    orders = []
    if C.degree() > 0:
        for f, e in C.factor():
            m = cyclotomic_order(f)
            if m == 0:
                raise RuntimeError("non-cyclotomic factor %s of cyclotomic part" % f)
            orders.extend([m] * e)
    orders.sort()
    if C[0] < 0:
        C = -C
    S = P // C
    return orders, C, S


def lcm_closure(ms):
    """
    The sorted list of least common multiples of subsets of ``ms``
    (including the empty subset, with lcm 1).
    """
    closure = {1}
    for m in sorted(set(ms)):
        closure |= {lcm(d, m) for d in closure}
    return sorted(closure)


def rank_over_extension(orders, n):
    """
    The number of roots of unity among the roots of the cyclotomic part
    that become 1 over the degree n extension: sum of phi(m) over the
    m dividing n.
    """
    return sum(euler_phi(m) for m in orders if n % m == 0)


def cyclotomic_artin_tate_factor(orders, n):
    """
    The product of (1 - zeta^n) over the roots zeta of the cyclotomic part
    with zeta^n != 1.
    """
    ans = ZZ(1)
    for m in orders:
        g = gcd(m, n)
        mp = m // g
        if mp == 1:
            continue
        # zeta^n runs over the primitive mp-th roots of unity,
        # each phi(m)/phi(mp) times
        ans *= ZZ(phi_at_one(mp)) ** (euler_phi(m) // euler_phi(mp))
    return ans
