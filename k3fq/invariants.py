"""
Invariants of a K3-type polynomial

    P(T) = q + a_1 T + ... + a_22 T^22,   a_22 = sign * q,

with all roots on the unit circle.  See schemas/k3_fq_isog.md for the
definitions of the columns computed here.
"""
from collections import Counter
from fractions import Fraction
from math import gcd

from sage.all import ZZ, ComplexField, RealField, euler_phi, moebius, pari, pi

from .cyclotomic import (R, T, split_cyclotomic, lcm, lcm_closure, rank_over_extension,
                         cyclotomic_artin_tate_factor)
from .galois import chebyshev_real_poly, galois_label
from .labels import poly_to_label, cremona_letter_code, DEGREE, HALF

NUM_COUNTS = 10  # number of point counts stored

# The 22-dimensional motive M sits inside the cohomology of each type of
# variety; outside M everything is spanned by Tate classes defined over F_q.
# For each type we record (weight shift w, [c_0, c_1, ...]) meaning
#   #Y(F_{q^n}) = sum_j c_j q^{jn} + q^{wn} Tr(Frob^n | M(-1)) / q^n
# i.e. the trace contribution is q^{(w-1)n} t_n with t_n = Tr(Frob^n | H^2).
VARIETY_TYPES = {
    # K3 surface: H^0 + H^2 + H^4, M = H^2(1)
    "k3": (1, [1, 0, 1]),
    # cubic fourfold: b_2 = 1, H^4 = Q h^2 + H^4_prim, M = H^4_prim(1)
    "cubic": (2, [1, 1, 1, 1, 1]),
    # Gushel-Mukai fourfold: b_2 = 1, H^4 = 2 Tate classes + H^4_van
    "gm": (2, [1, 1, 2, 1, 1]),
    # Verra fourfold: b_2 = 2, b_4 = 24 with three classes pulled back from
    # P^2 x P^2; we take M to contain one of them (to be confirmed)
    "verra": (2, [1, 2, 2, 2, 1]),
}


def newton_polygon(points):
    """
    Lower convex hull of ``points`` (a list of (x, y) with distinct
    increasing x, y Fractions); returns the list of vertices.
    """
    hull = []
    for pt in points:
        while len(hull) >= 2:
            (x1, y1), (x2, y2) = hull[-2], hull[-1]
            # remove hull[-1] if it lies on or above the segment hull[-2] -> pt
            if (y2 - y1) * (pt[0] - x1) >= (pt[1] - y1) * (x2 - x1):
                hull.pop()
            else:
                break
        hull.append(pt)
    return hull


def slopes_from_vertices(vertices):
    slopes = []
    for (x1, y1), (x2, y2) in zip(vertices, vertices[1:]):
        s = Fraction(y2 - y1) / (x2 - x1)
        slopes.extend([s] * (x2 - x1))
    return slopes


def polygon_value(vertices, x):
    """
    Value at ``x`` of the piecewise linear function with the given vertices.
    """
    for (x1, y1), (x2, y2) in zip(vertices, vertices[1:]):
        if x1 <= x <= x2:
            return y1 + Fraction(y2 - y1) * (x - x1) / (x2 - x1)
    raise ValueError("x out of range")


def ceil_frac(a):
    return -((-a.numerator) // a.denominator)


def slopes_text(slopes):
    """
    The abvar convention: each slope gets a letter A, B, ... according to
    how many times it has already appeared.
    """
    cntr = Counter()
    ans = []
    for s in slopes:
        ans.append(str(s) + cremona_letter_code(cntr[s]).upper())
        cntr[s] += 1
    return ans


def power_sums(Lcoeffs, nmax):
    """
    t_n = Tr(Frob^n) = sum of n-th powers of the reciprocal roots of
    L(T) = 1 + b_1 T + ..., via Newton's identities.
    """
    b = list(Lcoeffs) + [0] * max(0, nmax + 1 - len(Lcoeffs))
    t = [None]
    for n in range(1, nmax + 1):
        s = n * b[n]
        for i in range(1, n):
            s += b[i] * t[n - i]
        t.append(-s)
    return t


def closed_point_counts_nonnegative(N, nmax):
    """
    Whether sum_{d | n} mu(n/d) N[d] >= 0 for 1 <= n <= nmax.
    """
    for n in range(1, nmax + 1):
        s = sum(moebius(n // d) * N[d] for d in ZZ(n).divisors())
        if s < 0:
            return False
    return True


def _is_significant(rel, sprec):
    m = min(abs(c) for c in rel)
    M = max(abs(c) for c in rel)
    if (m + 1).exact_log(2) >= sprec:
        return False
    elif (M + 1).exact_log(2) >= sprec:
        raise RuntimeError("mixed significance in lindep")
    return True


def numerical_rank(numbers, sprec):
    """
    Dimension of the Q-span of ``numbers`` (real numbers with about sprec^2
    bits of precision), using integer relation finding as in abvar-fq.
    """
    numbers = list(numbers)
    while len(numbers) > 1:
        rel = [ZZ(c) for c in pari(numbers).lindep()]
        if not _is_significant(rel, sprec):
            return len(numbers)
        i = next(i for i, c in enumerate(rel) if c != 0)
        numbers.pop(i)
    return len(numbers)


class K3Poly:
    """
    A K3-type polynomial.  ``coeffs`` is the coefficient list (constant term
    first) of P(T), of length 23.
    """
    def __init__(self, coeffs):
        coeffs = [ZZ(c) for c in coeffs]
        if len(coeffs) != DEGREE + 1:
            raise ValueError("expected a polynomial of degree %s" % DEGREE)
        q = coeffs[0]
        if q < 2 or not q.is_prime_power():
            raise ValueError("constant term must be a prime power")
        sign = coeffs[-1] // q
        if sign not in (1, -1) or coeffs[-1] != sign * q:
            raise ValueError("leading coefficient must be +-q")
        if any(coeffs[DEGREE - i] != sign * coeffs[i] for i in range(HALF + 1)):
            raise ValueError("polynomial is not (anti)palindromic")
        self.coeffs = coeffs
        self.P = R(coeffs)
        self.q = q
        self.p, self.r = q.is_prime_power(get_data=True)
        self.sign = sign

    # ------------------------------------------------------------------
    # Factorization

    def _split(self):
        if not hasattr(self, "_orders"):
            self._orders, self._C, self._S = split_cyclotomic(self.P)
        return self._orders, self._C, self._S

    @property
    def cyclotomic_factors(self):
        return self._split()[0]

    @property
    def transcendental_poly(self):
        return self._split()[2]

    def transcendental_factors(self):
        """
        List of pairs (f, e), f irreducible primitive palindromic with
        f(0) > 0.
        """
        if hasattr(self, "_tfactors"):
            return self._tfactors
        S = self.transcendental_poly
        if S.degree() == 0:
            ans = []
        elif self.r == 1:
            # S is irreducible: a nontrivial factorization would have a
            # monic factor with roots on the unit circle, which is
            # cyclotomic by Kronecker's theorem.
            ans = [(S, 1)]
        else:
            ans = []
            for f, e in S.factor():
                if f.degree() == 0:
                    continue
                if f[0] < 0:
                    f = -f
                ans.append((f, e))
            ans.sort(key=lambda fe: (fe[0].degree(), fe[0].list()))
        self._tfactors = ans
        return ans

    # ------------------------------------------------------------------
    # Newton polygon

    def newton_data(self):
        """
        Vertices of the Newton polygon of L(T) = q^{-1} P(qT) in q-adic
        units: the points are (0, 0) and (i, i - 1 + v_p(a_i)/r).
        """
        p, r = self.p, self.r
        pts = [(0, Fraction(0))]
        for i in range(1, DEGREE + 1):
            a = self.coeffs[i]
            if a:
                pts.append((i, Fraction(i - 1) + Fraction(int(a.valuation(p)), int(r))))
        return newton_polygon(pts)

    # ------------------------------------------------------------------
    # Roots and angles

    def _transcendental_roots(self, prec):
        """
        For each distinct irreducible transcendental factor, its complex
        roots in the upper half plane at the given precision, with the
        multiplicity of the factor.
        """
        ans = []
        for f, e in self.transcendental_factors():
            roots = pari(f).polroots(precision=prec)
            CC = ComplexField(prec)
            up = [CC(z) for z in roots if CC(z).imag() > 0]
            assert len(up) == f.degree() // 2
            ans.append((up, e))
        return ans

    def angles(self):
        angles = []
        for m in self.cyclotomic_factors:
            if m == 1:
                angles.append(0.0)
            elif m == 2:
                angles.append(1.0)
            else:
                angles.extend(2.0 * j / m for j in range(1, (m + 1) // 2) if gcd(j, m) == 1)
        RRpi = RealField(53)(pi)
        for up, e in self._transcendental_roots(64):
            for z in up:
                angles.extend([float(z.arg() / RRpi)] * e)
        return sorted(angles)

    def angle_rank(self, certified):
        """
        ``certified`` says whether the transcendental part is f^e with f
        irreducible and Gal(f) the full wreath product C2 wr Gal(h); then
        the angle rank is deg(f)/2 (see galois.py).
        """
        tf = self.transcendental_factors()
        if not tf:
            return 0
        if certified and len(tf) == 1:
            return tf[0][0].degree() // 2
        sprec = 25
        prec = sprec ** 2
        RF = RealField(prec)
        angles = []
        for up, e in self._transcendental_roots(prec):
            angles.extend(z.arg() / RF(pi) for z in up)
        return numerical_rank(angles + [RF(1)], sprec) - 1

    # ------------------------------------------------------------------
    # Artin-Tate

    def artin_tate_value(self, n):
        """
        AT_n = (P_n(T) / (1-T)^{rho_n})(1), with P_n = q^n prod (1 - alpha^n T).
        """
        orders, C, S = self._split()
        cyc = cyclotomic_artin_tate_factor(orders, n)
        if S.degree() == 0:
            trans = self.q ** n
        else:
            # Res(S^rev, 1 - T^n) = q^n prod (1 - alpha_i^n), since S^rev has
            # leading coefficient S(0) = q and roots alpha_i.
            trans = S.reverse().resultant(1 - T ** n)
        return cyc * trans

    # ------------------------------------------------------------------

    def L_coeffs(self):
        """
        Coefficients of L(T) = 1 + sum q^{i-1} a_i T^i.
        """
        q = self.q
        return [ZZ(1)] + [q ** (i - 1) * self.coeffs[i] for i in range(1, DEGREE + 1)]

    def counts(self, kind, nmax):
        w, cs = VARIETY_TYPES[kind]
        q = self.q
        t = self._power_sums(nmax)
        return [None] + [sum(c * q ** (j * n) for j, c in enumerate(cs)) + q ** ((w - 1) * n) * t[n]
                         for n in range(1, nmax + 1)]

    def _power_sums(self, nmax):
        if getattr(self, "_t", None) is None or len(self._t) <= nmax:
            self._t = power_sums(self.L_coeffs(), nmax)
        return self._t

    def count_check_bound(self):
        """
        An n such that, for every variety type, closed point counts of degree
        > n are automatically positive: we use q^n >= 10^4 (see the README).
        """
        n = 1
        while self.q ** n < 10 ** 4:
            n += 1
        return n

    # ------------------------------------------------------------------

    def real_poly(self):
        if self.sign == 1:
            return chebyshev_real_poly(self.coeffs)
        Pprime = self.P // (1 - T ** 2)
        return chebyshev_real_poly(Pprime.list())

    def compute(self):
        """
        Return a dictionary of the columns of k3_fq_isog computed from the
        polynomial alone (no twist, base change or Kummer data).
        """
        q, p, r, sign = self.q, self.p, self.r, self.sign
        coeffs = self.coeffs
        orders, C, S = self._split()
        D = {}
        D["label"] = poly_to_label(coeffs)
        D["q"] = q
        D["p"] = p
        D["r"] = r
        D["sign"] = sign
        D["poly"] = list(coeffs)
        D["poly_str"] = " ".join(map(str, coeffs)) + " "
        D["real_poly"] = self.real_poly()

        # Newton polygon
        V = self.newton_data()
        slopes = slopes_from_vertices(V)
        D["slopes"] = slopes_text(slopes)
        h = sum(1 for s in slopes if s < 1)
        D["height"] = h if h else None
        if h == 0:
            k3np = True
        else:
            k3np = (slopes == [1 - Fraction(1, h)] * h + [Fraction(1)] * (DEGREE - 2 * h)
                    + [1 + Fraction(1, h)] * h)
        D["has_k3_newton_polygon"] = k3np
        elev = coelev = 0
        for x in range(0, HALF + 1):
            c = ceil_frac(polygon_value(V, x))
            if x >= 1:
                elev += max(0, c - (x - 1))
            coelev += max(0, x - c)
        D["newton_elevation"] = elev
        D["newton_coelevation"] = coelev
        D["is_ordinary"] = (h == 1)
        D["is_supersingular"] = (h == 0)

        # Cyclotomic part / Picard ranks / Artin-Tate
        rho = orders.count(1)
        rhobar = C.degree()
        assert rhobar == sum(euler_phi(m) for m in orders)
        assert (h == 0) == (rhobar == DEGREE)
        D["cyclotomic_factors"] = orders
        D["picard_rank"] = rho
        D["geometric_picard_rank"] = rhobar
        E = 1
        for m in orders:
            E = lcm(E, m)
        D["geometric_extension_degree"] = E
        ext = lcm_closure(orders)
        D["extension_degrees"] = ext
        D["picard_ranks"] = [rank_over_extension(orders, n) for n in ext]
        AT = [self.artin_tate_value(n) for n in ext]
        assert all(a > 0 for a in AT), (coeffs, AT)
        D["artin_tate_values"] = AT
        D["artin_tate_value"] = AT[0]
        D["ns_disc_class"] = (-1) ** (rho - 1) * AT[0].squarefree_part() if rho else None
        S1 = self.P // (1 - T) ** rho
        D["ej_condition"] = ZZ(S1(-1)).is_square()

        # Transcendental part
        D["transcendental_poly"] = S.list()
        tf = self.transcendental_factors()
        D["transcendental_factors"] = [[[int(c) for c in f.list()], int(e)] for f, e in tf]
        D["transcendental_irreducible"] = (len(tf) == 1 and tf[0][1] == 1) if tf else None
        gals = []
        certified = True
        for f, e in tf:
            lab, full = galois_label(f)
            gals.append(lab)
            certified = certified and full
        D["galois_groups"] = gals
        D["number_fields"] = None
        D["number_field_coeffs"] = None

        # Angles
        D["angles"] = self.angles()
        delta = self.angle_rank(certified)
        D["angle_rank"] = delta
        D["angle_corank"] = (DEGREE - rhobar) // 2 - delta

        # Point counts
        nmax = max(NUM_COUNTS, self.count_check_bound())
        growth = {}
        for kind in VARIETY_TYPES:
            N = self.counts(kind, nmax)
            D[kind + "_counts"] = N[1:NUM_COUNTS + 1]
            growth[kind] = closed_point_counts_nonnegative(N, nmax)
        N = self.counts("k3", nmax)
        D["k3_count"] = N[1]
        D["has_nonnegative_counts"] = all(c >= 0 for c in N[1:])
        D["has_count_growth"] = growth["k3"]
        D["is_pointless"] = (N[1] == 0)

        # Realizability: only necessary conditions at this stage
        k3_ok = rho >= 1 and D["ej_condition"] and growth["k3"] and k3np
        D["realizable_k3"] = 0 if k3_ok else -1
        D["realizable_nck3"] = 0
        D["realizable_cubic"] = 0 if (growth["cubic"] and k3np) else -1
        D["realizable_nccubic"] = 0
        D["realizable_gm"] = 0 if (growth["gm"] and k3np) else -1
        D["realizable_ncgm"] = 0
        D["realizable_verra"] = 0 if (growth["verra"] and k3np) else -1
        D["realizable_ncverra"] = 0
        return D
