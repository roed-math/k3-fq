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
from collections import defaultdict
from math import gcd

from sage.all import ZZ, QQ, RealField, ComplexField, PolynomialRing, pari, pi, euler_phi, lcm as sage_lcm

from .cyclotomic import R, T
from .labels import poly_to_label

M_TWIST = sage_lcm([m for m in range(1, 500) if euler_phi(m) <= 44])
PREC = 400
KEY_DIGITS = 30
MAX_TWISTS = 50


def base_change(coeffs, n):
    """
    Coefficients of P_n(T) = q^n prod (1 - alpha_i^n T), computed exactly as
    a resultant (as in abvar-fq's base_change).
    """
    P = R(coeffs)
    for ell, e in ZZ(n).factor():
        for _ in range(e):
            P = _base_change_prime(P, ell)
    return P.list()


_Ru = PolynomialRing(R, "u")


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
    All 22 root angles (fractions of a turn): exact rationals for the
    cyclotomic part, reals for the transcendental part.
    """
    cyc = []
    for m in row["cyclotomic_factors"]:
        cyc.extend(QQ((j, m)) for j in range(m) if gcd(j, m) == 1)
    return cyc, _angles(row["transcendental_factors"])


def _same_after(n, A, B, eps):
    """
    Whether the multisets {n a} and {n b} (mod 1) agree, for the exact parts
    and (within eps) for the numerical parts.
    """
    ca, ta = A
    cb, tb = B
    # (Sage's % on rationals is modular inversion, not the fractional part)
    if sorted(n * x - (n * x).floor() for x in ca) != sorted(n * x - (n * x).floor() for x in cb):
        return False

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


def _relevant_lcm(A, B, eps):
    """
    lcm of the orders of the roots of unity relating roots of A and B,
    together with the orders of the cyclotomic roots.
    """
    L = 1
    for x in A[0] + B[0]:
        L = sage_lcm(L, x.denominator())
    for a in A[1]:
        for b in B[1]:
            d = M_TWIST * (b - a)
            r = d.round()
            if abs(d - r) < eps * M_TWIST:
                L = sage_lcm(L, M_TWIST // gcd(int(r) % M_TWIST, M_TWIST))
    return ZZ(L)


# Base change labels are only recorded for q^n < 2^BC_LABEL_BITS
BC_LABEL_BITS = 128


def angle_data(row):
    """
    The data needed for twist computations: the geometric key and the
    angles of all roots.
    """
    A = _all_angles(row)
    return _key(row, A[1]), A


def minimal_twist_degrees(rowA, rowB, A=None, B=None):
    """
    The minimal n (for divisibility) with P_n = P'_n.  ``A`` and ``B`` are
    the outputs of ``_all_angles``.  Each n with q^n < 2^BC_LABEL_BITS is
    also checked by exact base change.
    """
    eps = RealField(PREC)(2) ** (-PREC // 2)
    if A is None:
        A = _all_angles(rowA)
    if B is None:
        B = _all_angles(rowB)
    L = _relevant_lcm(A, B, eps)
    bits = ZZ(rowA["q"]).nbits()
    found = []
    for n in L.divisors():
        if any(n % d == 0 for d in found):
            continue
        if _same_after(n, A, B, eps):
            if n * bits <= BC_LABEL_BITS and base_change(rowA["poly"], n) != base_change(rowB["poly"], n):
                raise RuntimeError("numerical twist check failed for %s, %s, n=%s"
                                   % (rowA["label"], rowB["label"], n))
            found.append(n)
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
