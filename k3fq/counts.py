"""
Counting K3-type polynomials by geometric Picard rank, without listing them.

A K3-type polynomial factors uniquely as C(T) S(T), with C a product of
cyclotomic polynomials (normalized to C(0) = 1) of degree 22 - 2k, and S a
transcendental part of degree 2k (constant and leading coefficient q, roots
on the unit circle, no cyclotomic factor).  The geometric Picard rank is
22 - 2k.  Hence

    N(q, rho >= R) = sum_{2k <= 22 - R} T(q, 2k) C(22 - 2k),

where C(e) counts products of cyclotomic polynomials of degree e (any
sign), and T(q, 2k) counts transcendental parts.

T comes from the counts W(q, 2k) printed by ``weilpoly 2k 1 q`` (all
palindromic polynomials with leading coefficient q and roots on the unit
circle, including those with cyclotomic factors).  Such a polynomial is
C S with C palindromic, i.e. with an even power of (1 - T), so

    W(q, 2k) = sum_{j <= k} T(q, 2j) C^+(2k - 2j),

where C^+ counts cyclotomic products with an even power of Phi_1.  This
can be solved for T recursively.
"""
import os
import re
import subprocess
import time

from sage.all import ZZ, PowerSeriesRing, euler_phi

from .cyclotomic import cyclotomic_orders_up_to
from .enumerate import WEILPOLY

DEGREE = 22


def cyclotomic_counts(dmax=DEGREE):
    """
    Lists C, Cplus with C[e] (resp. Cplus[e]) the number of products of
    cyclotomic polynomials of degree e (resp. with an even power of Phi_1).
    """
    S = PowerSeriesRing(ZZ, "x", default_prec=dmax + 1)
    x = S.gen()
    rest = S(1)
    for m in cyclotomic_orders_up_to(dmax):
        if m == 1:
            continue
        rest *= 1 / (1 - x ** euler_phi(m))
    C = (rest / (1 - x)).list()[:dmax + 1]
    Cplus = (rest / (1 - x ** 2)).list()[:dmax + 1]
    return [int(c) for c in C], [int(c) for c in Cplus]


def weil_count(d, q, threads=4):
    """
    The number of palindromic integer polynomials of degree d with leading
    and constant coefficient q and all roots on the unit circle, using
    Kedlaya's driver (output discarded).  Returns (count, seconds).
    """
    if d == 0:
        return 1, 0.0
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    t0 = time.time()
    res = subprocess.run([WEILPOLY, str(d), "1", str(q)], stdout=subprocess.DEVNULL,
                         stderr=subprocess.PIPE, env=env, text=True, check=True)
    m = re.search(r"Total number of polynomials: (\d+)", res.stderr)
    return int(m.group(1)), time.time() - t0


def transcendental_counts(W, Cplus):
    """
    Solve W[k] = sum_{j <= k} T[j] Cplus[2k - 2j] for T (indexed by k).
    """
    T = []
    for k, w in enumerate(W):
        T.append(w - sum(T[j] * Cplus[2 * k - 2 * j] for j in range(k)))
        assert T[-1] >= 0
    return T


def picard_table(T, C):
    """
    Given T[k] = T(q, 2k) for k <= K, return {R: N(q, rho >= R)} for the
    even R >= 22 - 2K.
    """
    out = {}
    total = 0
    for k, t in enumerate(T):
        total += t * C[DEGREE - 2 * k]
        out[DEGREE - 2 * k] = total
    return out


def parse_qs(specs):
    """
    Prime powers from arguments like "2", "3-64".
    """
    qs = []
    for spec in specs:
        if "-" in spec:
            a, b = map(int, spec.split("-"))
            qs.extend(q for q in range(a, b + 1) if ZZ(q).is_prime_power())
        else:
            qs.append(int(spec))
    return qs


def counts_table(qs, max_degree, threads=4, out=None):
    """
    Print a markdown table of N(q, rho >= R) for R >= 22 - max_degree.
    """
    C, Cplus = cyclotomic_counts()
    K = max_degree // 2
    Rs = [DEGREE - 2 * k for k in range(K + 1)]
    lines = ["| q | " + " | ".join("rho >= %s" % R for R in Rs) + " | seconds |",
             "| --- " * (len(Rs) + 2) + "|"]
    print("\n".join(lines), flush=True)
    results = {}
    for q in qs:
        W, secs = [], 0.0
        for k in range(K + 1):
            w, t = weil_count(2 * k, q, threads)
            W.append(w)
            secs += t
        tab = picard_table(transcendental_counts(W, Cplus), C)
        results[q] = tab
        line = "| %s | %s | %.1f |" % (q, " | ".join(str(tab[R]) for R in Rs), secs)
        lines.append(line)
        print(line, flush=True)
    if out:
        with open(out, "w") as F:
            F.write("\n".join(lines) + "\n")
    return results
