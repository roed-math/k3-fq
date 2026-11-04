"""
Twists and base change.

Two K3-type polynomials P, P' over F_q are twists if P_n = P'_n for some n,
where P_n(T) = q^n prod (1 - alpha_i^n T) is the base change to F_{q^n}.

Geometric key.  If P_n = P'_n then each transcendental root beta of P' is
zeta * alpha for a root alpha of P and a root of unity zeta.  Here zeta is an
eigenvalue of a finite-order element of a Q-algebra acting on the 22-dimensional
realization, so phi(ord zeta) <= 22.  We take
M = lcm{m : phi(m) <= 44} (about 1.3e20) for safety.  Then the multiset of
{M theta} (theta = arg(alpha)/2pi, alpha over the transcendental roots,
fractional part) is invariant under twisting, and it determines the base
change to F_{q^M}, hence the twist class.  We compute it numerically at 400
bits and round to 30 digits.

Within a twist class, the set of n with P_n = P'_n is closed under multiples
and under gcd with L, the lcm of the orders of the roots of unity involved.
So its minimal elements are divisors of L.  We find them numerically and then
check each one by exact base change.
"""
from collections import Counter, defaultdict
from functools import lru_cache
from math import gcd

from sage.all import ZZ, QQ, RealField, ComplexField, PolynomialRing, pari, pi, euler_phi, lcm as sage_lcm

from .cyclotomic import R, T
from .labels import poly_to_label
from .invariants import power_sums

M_TWIST = sage_lcm([m for m in range(1, 500) if euler_phi(m) <= 44])
PREC = 400
KEY_DIGITS = 30
MAX_TWISTS = 50


def base_change(coeffs, n):
    """
    Coefficients of P_n(T) = q^n prod (1 - alpha_i^n T).

    Write L(T) = q^{-1} P(qT), whose reciprocal roots qalpha_i have power sums
    t_j = Tr(Frob^j).  The base change L_n has power sums t_{kn}, so Newton's
    identities give L_n, and P_n(T) = q^n L_n(T / q^n).  Everything is exact
    integer arithmetic.
    """
    n = int(n)
    q = ZZ(coeffs[0])
    d = len(coeffs) - 1
    if n == 1:
        return [ZZ(c) for c in coeffs]
    b = [ZZ(1)] + [q ** (i - 1) * ZZ(coeffs[i]) for i in range(1, d + 1)]
    t = power_sums(b, d * n)
    u = [None] + [t[k * n] for k in range(1, d + 1)]
    c = [ZZ(1)]
    for k in range(1, d + 1):
        s = sum(u[i] * c[k - i] for i in range(1, k + 1))
        assert s % k == 0
        c.append(-s // k)
    qn = q ** n
    ans = [qn]
    for i in range(1, d + 1):
        num = c[i]
        den = qn ** (i - 1)
        assert num % den == 0
        ans.append(num // den)
    return ans


_Ru = PolynomialRing(R, "u")


def base_change_resultant(coeffs, n):
    """
    The same as ``base_change``, via resultants (slower; used in tests).
    """
    P = R(coeffs)
    for ell, e in ZZ(n).factor():
        for _ in range(e):
            P = _base_change_prime(P, ell)
    return P.list()


def _base_change_prime(P, ell):
    """
    One step of base change, of prime degree ell.  The resultant
    Res_u(P(u), u^ell - T) has roots T = beta^ell for the roots beta of P;
    we normalize it to have constant term P(0)^ell.
    """
    u = _Ru.gen()
    res = P(u).resultant(u ** ell - T)
    q = ZZ(P[0])
    c0 = ZZ(res[0])
    Pn = [ZZ(c) * q ** ell for c in res.list()]
    assert all(c % c0 == 0 for c in Pn)
    return R([c // c0 for c in Pn])


def _angles(transcendental_factors, prec=PREC):
    """
    The arguments (as fractions of a full turn, in [0,1)) of all roots of
    the transcendental part, with multiplicity.
    """
    RF = RealField(prec)
    CF = ComplexField(prec)
    twopi = 2 * RF(pi)
    ans = []
    for f, e in transcendental_factors:
        for z in pari(R(f)).polroots(precision=prec):
            theta = CF(z).arg() / twopi
            ans.extend([theta - theta.floor()] * e)
    return ans


def geometric_key(row):
    """
    A string determining the twist class of ``row`` (a dictionary with at
    least q, geometric_picard_rank and transcendental_factors).
    """
    return _key(row, _angles(row["transcendental_factors"]))


def _key(row, thetas):
    RF = RealField(PREC)
    vals = []
    scale = RF(10) ** KEY_DIGITS
    for theta in thetas:
        x = M_TWIST * theta
        x = x - x.floor()
        v = (x * scale).round()
        if v == scale:
            v = 0
        vals.append(int(v))
    vals.sort()
    return "%s:%s:%s" % (row["q"], row["geometric_picard_rank"], ",".join(map(str, vals)))


def _all_angles(row):
    """
    The data used to compare base changes: the cyclotomic orders (as a
    Counter) and the angles of the transcendental roots.
    """
    return Counter(row["cyclotomic_factors"]), _angles(row["transcendental_factors"])


@lru_cache(maxsize=None)
def _phi(m):
    return int(euler_phi(m))


def _cyclotomic_base_change(orders, n):
    """
    The cyclotomic orders of the base change to degree n: a root of unity of
    order m becomes one of order m/gcd(m, n), and the phi(m) roots of
    Phi_m map onto the primitive (m/g)-th roots phi(m)/phi(m/g) times each.
    """
    out = Counter()
    for m, e in orders.items():
        mp = m // gcd(m, n)
        out[mp] += e * (_phi(m) // _phi(mp))
    return out


def _transcendental_same_after(n, ta, tb, eps):
    """
    Whether the multisets {n a} and {n b} (mod 1) agree within eps.
    """
    def red(xs):
        out = []
        for x in xs:
            y = n * x
            y = y - y.floor()
            if 1 - y < eps:
                y = y - 1
            out.append(y)
        return sorted(out)
    return all(abs(x - y) < eps for x, y in zip(red(ta), red(tb)))


def _transcendental_lcm(ta, tb, eps):
    """
    lcm of the orders of the roots of unity b/a, over pairs of transcendental
    roots whose ratio is a root of unity.
    """
    L = 1
    for a in ta:
        for b in tb:
            d = M_TWIST * (b - a)
            r = d.round()
            if abs(d - r) < eps * M_TWIST:
                L = sage_lcm(L, M_TWIST // gcd(int(r) % M_TWIST, M_TWIST))
    return ZZ(L)


def _minimal_divisors(L, valid):
    """
    The minimal divisors n of L (for divisibility) with valid(n), assuming
    the set of valid n is closed under taking multiples.
    """
    found = []
    for n in ZZ(L).divisors():
        if any(n % d == 0 for d in found):
            continue
        if valid(n):
            found.append(n)
    return found


def _minimal_elements(ns):
    ns = sorted(set(ns))
    return [n for i, n in enumerate(ns) if not any(n % d == 0 for d in ns[:i])]


# Base change labels are only recorded for q^n < 2^BC_LABEL_BITS
BC_LABEL_BITS = 128


def angle_data(row):
    """
    The data needed for twist computations: the geometric key and the
    output of ``_all_angles``.
    """
    A = _all_angles(row)
    return _key(row, A[1]), A


def minimal_twist_degrees(rowA, rowB, A=None, B=None):
    """
    The minimal n (for divisibility) with P_n = P'_n.  ``A`` and ``B`` are
    the outputs of ``_all_angles``.

    The n for which the cyclotomic parts agree, and those for which the
    transcendental parts agree, are both closed under multiples.  So the
    minimal n for both are minimal elements of lcms of the minimal ones.
    Their minimal elements divide the lcm of the cyclotomic orders, and of the
    orders of the ratios of transcendental roots, respectively.  Each n with
    q^n < 2^BC_LABEL_BITS is also checked by exact base change.
    """
    eps = RealField(PREC)(2) ** (-PREC // 2)
    if A is None:
        A = _all_angles(rowA)
    if B is None:
        B = _all_angles(rowB)
    (ca, ta), (cb, tb) = A, B
    LC = 1
    for m in list(ca) + list(cb):
        LC = sage_lcm(LC, m)
    cyc = _minimal_divisors(LC, lambda n: _cyclotomic_base_change(ca, n) == _cyclotomic_base_change(cb, n))
    trans = _minimal_divisors(_transcendental_lcm(ta, tb, eps),
                              lambda n: _transcendental_same_after(n, ta, tb, eps))
    found = _minimal_elements(sage_lcm(c, t) for c in cyc for t in trans)
    bits = ZZ(rowA["q"]).nbits()
    for n in found:
        if n * bits <= BC_LABEL_BITS and base_change(rowA["poly"], n) != base_change(rowB["poly"], n):
            raise RuntimeError("numerical twist check failed for %s, %s, n=%s"
                               % (rowA["label"], rowB["label"], n))
    return found


def class_twists(members, angles):
    """
    ``members``: rows in one twist class (at most MAX_TWISTS of them);
    ``angles``: label -> output of ``_all_angles``.  Returns label -> list of
    [label', bc_label, n], with bc_label None when q^n >= 2^BC_LABEL_BITS.
    """
    twists = {r["label"]: [] for r in members}
    bits = ZZ(members[0]["q"]).nbits()
    for i, A in enumerate(members):
        for B in members[i + 1:]:
            for d in minimal_twist_degrees(A, B, angles[A["label"]], angles[B["label"]]):
                d = int(d)
                bc = poly_to_label(base_change(A["poly"], d)) if d * bits <= BC_LABEL_BITS else None
                twists[A["label"]].append([B["label"], bc, d])
                twists[B["label"]].append([A["label"], bc, d])
    for lab in twists:
        twists[lab].sort(key=lambda t: (t[2], t[0]))
    return twists


def twist_columns(members, twists=None):
    """
    The TWIST_COLUMNS for each member of a twist class.
    """
    rep = min(r["label"] for r in members)
    out = {}
    for r in members:
        lab = r["label"]
        tw = None if twists is None else twists[lab]
        out[lab] = {
            "geometric_class": rep,
            "twist_count": len(members),
            "twists": tw,
            "max_twist_degree": None if tw is None else max((t[2] for t in tw), default=1),
        }
    return out


def twist_data(rows):
    """
    Given all rows over one F_q (dictionaries with label, q, poly,
    cyclotomic_factors, geometric_picard_rank, transcendental_factors),
    return a dictionary label -> dict of TWIST_COLUMNS.  (In memory; the
    twists stage does the same thing chunk by chunk.)
    """
    classes = defaultdict(list)
    angles = {}
    for row in rows:
        key, angles[row["label"]] = angle_data(row)
        classes[key].append(row)
    out = {}
    for members in classes.values():
        tw = class_twists(members, angles) if len(members) <= MAX_TWISTS else None
        out.update(twist_columns(members, tw))
    return out


def primitive_model_data(rows_by_q):
    """
    ``rows_by_q`` maps q to the list of rows (with label and poly) over F_q.
    Returns label -> {is_primitive, primitive_models} for every row, with
    is_primitive None when some proper subfield has no data.
    """
    out = {}
    for q, rows in rows_by_q.items():
        p, r = ZZ(q).is_prime_power(get_data=True)
        models = defaultdict(list)
        complete = True
        for d in ZZ(r).divisors():
            if d == 1:
                continue
            q0 = p ** (r // d)
            if q0 not in rows_by_q:
                complete = False
                continue
            for row0 in rows_by_q[q0]:
                models[poly_to_label(base_change(row0["poly"], d))].append(row0["label"])
        for row in rows:
            ms = sorted(models.get(row["label"], []))
            if ms:
                out[row["label"]] = {"is_primitive": False, "primitive_models": ms}
            else:
                out[row["label"]] = {"is_primitive": True if complete else None,
                                     "primitive_models": None}
    return out
