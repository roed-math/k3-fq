"""
Enumeration of K3-type polynomials using Kedlaya's root-unitary code.

We call the OpenMP driver ``c-script/weilpoly.c`` from the root-unitary
repository (built as ``bin/weilpoly`` by the Makefile).  ``weilpoly d 1 q``
prints all integer polynomials of degree d with leading and constant
coefficient q, sign +1 and all roots on the unit circle, one per line.

* Sign +1: ``weilpoly 22 1 q``.
* Sign -1: such a polynomial is (1 - T^2) P' with P' of degree 20 and sign
  +1 (this is how ``WeilPolynomials(..., sign=-1)`` works internally).
* Elevated geometric Picard rank >= R: as in the Kedlaya-Sutherland k3
  scripts, enumerate the transcendental parts of degree 2k <= 22 - R (the
  output of ``weilpoly 2k 1 q`` without cyclotomic factors) and multiply by
  all products of cyclotomic polynomials of degree 22 - 2k.
"""
import os
import subprocess
import sys
import time

from sage.all import ZZ, euler_phi, cyclotomic_polynomial

from .cyclotomic import R, T, cyclotomic_orders_up_to

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WEILPOLY = os.path.join(ROOT, "bin", "weilpoly")
DEGREE = 22


def run_weilpoly(d, q, threads):
    """
    Iterate over the coefficient lists printed by ``weilpoly d 1 q``.
    """
    if not os.path.exists(WEILPOLY):
        raise RuntimeError("%s not found: run make first" % WEILPOLY)
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    proc = subprocess.Popen([WEILPOLY, str(d), "1", str(q)], stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, env=env, text=True, bufsize=1 << 20)
    for line in proc.stdout:
        yield [int(c) for c in line.split()]
    err = proc.stderr.read()
    if proc.wait() != 0:
        raise RuntimeError("weilpoly failed: %s" % err)
    sys.stderr.write(err)


def _write_lines(filename, polys):
    tmp = filename + ".tmp"
    count = 0
    with open(tmp, "w") as F:
        for c in polys:
            F.write(" ".join(map(str, c)) + "\n")
            count += 1
    os.replace(tmp, filename)
    return count


def times_one_minus_T2(c):
    """
    Coefficients of (1 - T^2) * sum c_i T^i.
    """
    c = list(c) + [0, 0]
    return [c[i] - (c[i - 2] if i >= 2 else 0) for i in range(len(c))]


def enumerate_full(q, outdir, threads=8):
    """
    Write all K3-type polynomials over F_q to ``outdir/q{q}_p.txt`` (sign +1)
    and ``outdir/q{q}_m.txt`` (sign -1).  Returns the two counts.
    """
    os.makedirs(outdir, exist_ok=True)
    counts = []
    for s, gen in (("p", run_weilpoly(DEGREE, q, threads)),
                   ("m", (times_one_minus_T2(c) for c in run_weilpoly(DEGREE - 2, q, threads)))):
        t0 = time.time()
        n = _write_lines(os.path.join(outdir, "q%s_%s.txt" % (q, s)), gen)
        print("q=%s sign %s: %s polynomials (%.1fs)" % (q, s, n, time.time() - t0))
        counts.append(n)
    return counts


def cyclotomic_products(d):
    """
    All products of cyclotomic polynomials of total degree d, each
    normalized to have constant term 1, as polynomials in T.
    """
    ms = [m for m in cyclotomic_orders_up_to(max(d, 1)) if euler_phi(m) <= d]
    phis = {m: (1 - T if m == 1 else R(cyclotomic_polynomial(m).list())) for m in ms}
    degs = {m: euler_phi(m) for m in ms}

    def rec(i, remaining):
        if remaining == 0:
            yield R(1)
            return
        if i == len(ms):
            return
        m = ms[i]
        f, e = phis[m], degs[m]
        power = R(1)
        for k in range(remaining // e + 1):
            for rest in rec(i + 1, remaining - k * e):
                yield power * rest
            power *= f
    return list(rec(0, d))


def transcendental_parts(two_k, q, threads):
    """
    Polynomials of degree 2k with leading and constant coefficient q, roots
    on the unit circle and no cyclotomic factor.
    """
    if two_k == 0:
        return [R([q])]
    return [f for f in (R(c) for c in run_weilpoly(two_k, q, threads))
            if not f.has_cyclotomic_factor()]


def enumerate_elevated(q, min_rank, outdir, threads=8):
    """
    Write all K3-type polynomials over F_q with geometric Picard rank at
    least ``min_rank`` to ``outdir/q{q}_rho{min_rank}.txt``.
    """
    os.makedirs(outdir, exist_ok=True)

    def gen():
        for two_k in range(0, DEGREE - min_rank + 1, 2):
            t0 = time.time()
            S = transcendental_parts(two_k, q, threads)
            C = cyclotomic_products(DEGREE - two_k)
            print("q=%s degree %s: %s transcendental parts, %s cyclotomic parts (%.1fs)"
                  % (q, two_k, len(S), len(C), time.time() - t0))
            for s in S:
                for c in C:
                    P = c * s
                    assert P.degree() == DEGREE and P[0] == q
                    yield P.list()
    filename = os.path.join(outdir, "q%s_rho%s.txt" % (q, min_rank))
    n = _write_lines(filename, gen())
    print("q=%s, geometric Picard rank >= %s: %s polynomials" % (q, min_rank, n))
    return n
