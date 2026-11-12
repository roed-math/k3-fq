"""
The refined Artin-Tate conditions of Auel and Shaikh, "A cyclotomic Artin-Tate
formula and sphere packing for K3 surfaces" (arXiv:2611.13468, Theorem 4.1),
with the discriminant bounds of Auel, Engel, Kedlaya and Shaikh, "Pointless K3
surfaces over finite fields" (arXiv:2611.08356, Proposition 3.1).

For a K3 surface X over F_q and m >= 2, write

    P_X(t) = (1 - t)^rho C_m(t) P_m(t),   C_m = prod_{d | m, d > 1} Phi_d^{a_d},

with P_m(zeta_d) != 0 for d | m.  Then

    Res((t^m - 1)/(t - 1), P_m) = |Br_m(X)| |NS_m(X)| / |A_m(X)| > 0,

where NS_m(X), the Frobenius coinvariant lattice of NS(X_m), is (up to sign)
even, positive definite and rootless of rank deg C_m, and |A_m(X)| divides
gcd(m^rho, C_m(1)).  So when deg C_m > 0,

    Res > 0   and   gcd(m^rho, C_m(1)) Res >= delta_{deg C_m},

where delta_k is a lower bound for the discriminant of an even rootless
positive definite lattice of rank k.  The conditions also apply to every base
change X_l, whose Weil polynomial is that of X with its roots raised to the
l-th power.  We check the base changes l | E with l < E (E is the degree over
which all of NS is defined; for l = E there is nothing to check) and, for X_l,
every m | E/l with m >= 2.  (For other m, the conditions follow from those for
gcd(m, E/l).)
"""
from functools import lru_cache

from sage.all import ZZ, PolynomialRing, cyclotomic_polynomial, divisors, euler_phi, gcd

from .cyclotomic import phi_at_one
from .twists import base_change

R = PolynomialRing(ZZ, "t")
t = R.gen()

# Recorded in the pipeline's markers: change it when the conditions or bounds change.
VERSION = "Auel-Shaikh Thm 4.1, l | E, m | E/l; delta from Auel-Engel-Kedlaya-Shaikh Table 2"

# Lower bounds for the discriminant of an even rootless positive definite lattice of rank k
# (Auel-Engel-Kedlaya-Shaikh, Table 2; optimal for k <= 9 and k = 16).
DELTA = {1: 4, 2: 12, 3: 32, 4: 64, 5: 128, 6: 192, 7: 256, 8: 256, 9: 512, 10: 370, 11: 400,
         12: 410, 13: 406, 14: 370, 15: 318, 16: 256, 17: 194, 18: 134, 19: 86, 20: 50,
         21: 26, 22: 12, 23: 4}


@lru_cache(maxsize=None)
def _f(m):
    return (t ** m - 1) // (t - 1)


@lru_cache(maxsize=None)
def _cyclotomic_resultant(d, m):
    """Res((t^m - 1)/(t - 1), Phi_d) for d not dividing m (a nonzero integer)."""
    return _f(m).resultant(R(cyclotomic_polynomial(d).list()))


@lru_cache(maxsize=2 ** 16)
def _trans_base_change(scoeffs, l):
    if len(scoeffs) == 1:            # geometric Picard rank 22: S = q
        return (scoeffs[0] ** l,)
    return tuple(int(c) for c in base_change(list(scoeffs), l))


@lru_cache(maxsize=2 ** 18)
def _trans_resultant(scoeffs, m):
    return _f(m).resultant(R(list(scoeffs)))


def base_change_orders(orders, l):
    """The orders of the cyclotomic factors after base change to degree l."""
    out = []
    for d in orders:
        dl = d // gcd(d, l)
        out.extend([dl] * (euler_phi(d) // euler_phi(dl)))
    return sorted(out)


def refined_artin_tate_failures(orders, scoeffs, first=True):
    """
    The pairs (l, m) for which the class with cyclotomic factors ``orders``
    (with multiplicity, 1 included) and transcendental part ``scoeffs``
    (constant term q first) fails the refined Artin-Tate condition; with
    ``first``, stop at the first one.  Empty for classes that pass.
    """
    orders = [int(d) for d in orders]
    E = 1
    for d in orders:
        E = E * d // gcd(E, d)
    scoeffs = tuple(int(c) for c in scoeffs)
    failures = []
    for l in divisors(E):
        if l == E:
            continue
        ol = base_change_orders(orders, l) if l > 1 else sorted(orders)
        rho = ol.count(1)
        sl = _trans_base_change(scoeffs, l) if l > 1 else scoeffs
        for m in divisors(E // l):
            if m < 2:
                continue
            inside = [d for d in ol if d > 1 and m % d == 0]
            if not inside:
                continue
            k = sum(euler_phi(d) for d in inside)
            res = _trans_resultant(sl, m)
            for d in ol:
                if d > 1 and m % d != 0:
                    res *= _cyclotomic_resultant(d, m)
            c1 = ZZ(1)
            for d in inside:
                c1 *= phi_at_one(d)
            if res <= 0 or gcd(ZZ(m) ** rho, c1) * res < DELTA[k]:
                failures.append((int(l), int(m)))
                if first:
                    return failures
    return failures
