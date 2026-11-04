"""
Galois groups of the transcendental factors.

An irreducible transcendental factor f of a K3-type polynomial has degree 2k
and is palindromic with no real roots, so f(T) = T^k h(T + 1/T) for a real
polynomial h of degree k, and Gal(f) is a subgroup of the wreath product
C2 wr Gal(h) acting on 2k points with blocks {alpha, 1/alpha}.  Let N be the
kernel of Gal(f) -> Gal(h); it is an F_2[Gal(h)]-submodule of F_2^k.

Certificate.  Suppose f mod ell is squarefree of full degree and factors with
exactly one factor of degree 2 and all others of odd degree.  The Frobenius
sigma then has a fixed block that it swaps and odd-length cycles that do not
swap.  Raising sigma to the odd lcm of its cycle lengths gives the element of
N with exactly one nonzero coordinate.  Gal(h) is transitive, so N = F_2^k,
and Gal(f) is the full wreath product C2 wr Gal(h), with |Gal(f)| = 2^k |Gal(h)|.

In that case the reflection-type part of the permutation representation is
irreducible, and the angle rank is k (see invariants.py).

Gal(h) (degree <= 11) is computed with PARI's polgalois, which returns
transitive group numbers in the standard (GAP/Magma/LMFDB) numbering.  The
transitive label of C2 wr kTj is read from data/wreath_labels.txt, which
magma/wreath_labels.m generates.
"""
import os
from collections import Counter
from functools import lru_cache
from math import gcd

from sage.all import ZZ, pari, prime_range

_DATA = os.path.join(os.path.dirname(__file__), "data")

# polgalois only uses the standard transitive group numbering in the new format
pari.default("new_galois_format", 1)


@lru_cache(maxsize=None)
def wreath_labels():
    table = {}
    with open(os.path.join(_DATA, "wreath_labels.txt")) as F:
        for line in F:
            k, j, n, i = map(int, line.split())
            table[k, j] = (n, i)
    return table


def chebyshev_real_poly(coeffs):
    """
    Given the coefficient list (constant first) of a palindromic polynomial
    f of degree 2k (f_i = f_{2k-i}), return the coefficient list of h with
    f(T) = T^k h(T + 1/T).

    Uses T^j + T^-j = D_j(T + 1/T) with D_0 = 2, D_1 = x,
    D_j = x D_{j-1} - D_{j-2}.
    """
    d = len(coeffs) - 1
    if d % 2:
        raise ValueError("degree must be even")
    k = d // 2
    h = [0] * (k + 1)
    h[0] = coeffs[k]
    Dprev, D = [2], [0, 1]  # D_0, D_1
    for j in range(1, k + 1):
        c = coeffs[k + j]
        if c:
            for i, dc in enumerate(D):
                h[i] += c * dc
        # D_{j+1} = x D_j - D_{j-1}
        Dnext = [0] + D
        for i, dc in enumerate(Dprev):
            Dnext[i] -= dc
        Dprev, D = D, Dnext
    return h


def _degrees_mod(fpari, ell):
    """
    Degrees of the irreducible factors of f mod ell (f assumed squarefree
    mod ell), via distinct degree factorization.
    """
    M = pari.factormodDDF(fpari, ell)
    degs = []
    for i in range(len(M[0])):
        g, d = M[0][i], int(M[1][i])
        degs.extend([d] * (int(pari.poldegree(g)) // d))
    return degs


def frobenius_kernel_weight(hdegs, fdegs):
    """
    Let sigma be a Frobenius element, acting on the k blocks {alpha, 1/alpha}
    by a permutation pi with cycle lengths ``hdegs`` (the factor degrees of
    h mod ell).  A cycle of pi of length c either does not swap (two factors
    of f of degree c) or swaps an odd number of times (one factor of degree
    2c).  Let m be the lcm of the cycle lengths.  Then sigma^m lies in N and
    is supported on the cycles of length c that swap with m/c odd; return
    its weight.
    """
    a = Counter(hdegs)
    b = Counter(fdegs)
    odd = {}
    for c in sorted(a):
        # b_c = 2 * (non-swapping cycles of length c) + (swapping cycles of length c/2)
        from_half = odd.get(c // 2, 0) if c % 2 == 0 else 0
        even_c = (b[c] - from_half) // 2
        odd[c] = a[c] - even_c
    m = 1
    for c in a:
        m = m * c // gcd(m, c)
    return sum(c * odd[c] for c in a if (m // c) % 2 == 1)


def certify_full_wreath(f, h, H, max_primes=1000):
    """
    Look for a prime giving a certificate that N = F_2^k (module docstring).
    ``H`` = (k, j) identifies Gal(h) as kTj.  Returns the prime, or None.

    For any transitive Gal(h), an element of N of weight 1 suffices.  If
    Gal(h) is S_k or A_k (k >= 3), the only submodules of F_2^k are 0,
    <(1,...,1)>, the sum-zero submodule E and F_2^k; N contains (1,...,1)
    (complex conjugation).  So any element of N other than 0 and (1,...,1)
    gives N >= E, hence N = F_2^k when k is odd, or when the element has
    odd weight.
    """
    k = f.degree() // 2
    if k == 0:
        return None
    big = k >= 3 and H[1] in symmetric_or_alternating(k)
    lead = ZZ(f.leading_coefficient())
    disc = ZZ(f.discriminant())
    fpari = pari(f)
    hpari = pari(h)
    count = 0
    for ell in prime_range(3, 10**6):
        if lead % ell == 0 or disc % ell == 0:
            continue
        w = frobenius_kernel_weight(_degrees_mod(hpari, ell), _degrees_mod(fpari, ell))
        if w == 1 or (big and 0 < w < k and (k % 2 == 1 or w % 2 == 1)):
            return ell
        count += 1
        if count >= max_primes:
            return None
    return None


@lru_cache(maxsize=None)
def symmetric_or_alternating(k):
    """
    The transitive group numbers j such that kTj is S_k or A_k: the two
    largest groups in the standard numbering.
    """
    n = _NUM_TRANSITIVE[k]
    return (n - 1, n)


# number of transitive groups of degree k (k <= 11)
_NUM_TRANSITIVE = [0, 1, 1, 2, 5, 5, 16, 7, 50, 34, 45, 8]


def real_galois_group(h):
    """
    The pair (k, j) such that Gal(h) is kTj (h irreducible of degree k <= 11).
    """
    k = h.degree()
    if k == 1:
        return (1, 1)
    G = pari(h).polgalois()
    return (k, int(G[2]))


def galois_label(f, real_coeffs=None):
    """
    Try to determine the Galois group of the irreducible palindromic
    factor ``f`` (degree 2k, no real roots) as a transitive group label.

    Returns ``(label, full_wreath)`` where ``full_wreath`` is True if the
    certificate succeeded (and then ``label`` is that of C2 wr Gal(h)),
    and ``(None, False)`` if it failed (use magma/galois_groups.m).
    """
    k = f.degree() // 2
    if real_coeffs is None:
        real_coeffs = chebyshev_real_poly(f.list())
    h = f.parent()(real_coeffs)
    if 2 * k <= 11:
        # small enough for polgalois directly
        G = pari(f).polgalois()
        Horder = 1 if k == 1 else int(pari(h).polgalois()[0])
        return "%sT%s" % (2 * k, int(G[2])), int(G[0]) == 2**k * Horder
    H = real_galois_group(h)
    if certify_full_wreath(f, h, H) is None:
        return None, False
    n, i = wreath_labels()[H]
    assert n == 2 * k
    return "%sT%s" % (n, i), True
