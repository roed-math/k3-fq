"""
Input data for the lattice enumeration: one JSON line per class over F_q with
geometric Picard rank RHOBAR and Picard rank >= 1.

    sage -python lattices/prep.py TRANSFILE Q RHOBAR NSAMPLE SEED OUT

TRANSFILE lists transcendental parts of degree 22 - RHOBAR (one coefficient
list per line); NSAMPLE of them are taken at random (0 for all), times every
cyclotomic part of degree RHOBAR containing (1 - T).  Over F_2 the classes of
the quartic census with this geometric rank are always included.

Per class: label, rho, the orders d > 1 of the cyclotomic factors (with
multiplicity), E, realizable_k3 from the basic necessary conditions, and for
each m | E the Artin-Tate value AT_m = P_{X_m}^*(1) (= |Br X_m| |disc NS(X_m)|)
and, for m > 1, R_m = Res((t^m - 1)/(t - 1), P_m), where
P_X = (1 - t)^rho C_m(t) P_m(t) and C_m collects the Phi_d with d | m, d > 1
(Auel-Shaikh, Theorem 4.1).
"""
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from sage.all import ZZ, PolynomialRing, cyclotomic_polynomial, divisors

from k3fq.enumerate import cyclotomic_products
from k3fq.invariants import K3Poly
from k3fq.labels import poly_to_label, label_to_poly
from k3fq.stages import quartic_census_labels

R = PolynomialRing(ZZ, "t")
t = R.gen()


def class_data(coeffs, census=False):
    P = R(coeffs)
    D = K3Poly(coeffs).compute()
    rho = int(D["picard_rank"])
    orders = [int(d) for d in D["cyclotomic_factors"] if d != 1]
    E = int(D["geometric_extension_degree"])
    K = K3Poly(coeffs)
    AT = {}
    Res = {}
    for m in divisors(E):
        AT[int(m)] = str(K.artin_tate_value(m))
        if m > 1:
            Cm = R(1)
            for d in orders:
                if m % d == 0:
                    Cm *= cyclotomic_polynomial(d, t)
            Pm = P // ((1 - t) ** rho * Cm)
            assert Pm * (1 - t) ** rho * Cm == P
            Res[int(m)] = str(abs(((t ** m - 1) // (t - 1)).resultant(Pm)))
    return {"label": D["label"], "q": int(D["q"]), "rhobar": int(D["geometric_picard_rank"]), "rho": rho,
            "orders": orders, "E": E, "realizable_k3": int(D["realizable_k3"]),
            "ordinary": bool(D["is_ordinary"]), "census": census, "AT": AT, "Res": Res,
            "poly": [int(c) for c in coeffs], "S": [int(c) for c in D["transcendental_poly"]]}


def main(transfile, q, rhobar, nsample, seed, out):
    with open(transfile) as F:
        trans = [[int(c) for c in line.split()] for line in F]
    if nsample:
        trans = random.Random(seed).sample(trans, nsample)
    cycs = [c for c in cyclotomic_products(rhobar) if c(1) == 0]
    assert len(cycs) == 10 or rhobar != 4
    census = set()
    if q == 2:
        census = {lab for lab in quartic_census_labels() if lab.split(".")[1] == str(rhobar)}
    seen = set()
    with open(out, "w") as F:
        def emit(coeffs):
            lab = poly_to_label(coeffs)
            if lab in seen:
                return
            seen.add(lab)
            F.write(json.dumps(class_data(coeffs, lab in census)) + "\n")
        for s in trans:
            for c in cycs:
                P = R(c.list()) * R(s)
                assert P[0] == q
                emit([int(a) for a in P.list()])
        for lab in sorted(census):
            coeffs = label_to_poly(lab)
            if R(coeffs)(1) == 0:      # Picard rank >= 1 (always, for a K3 surface)
                emit(coeffs)
    print("%s classes (%s cyclotomic parts, %s census classes)" % (len(seen), len(cycs), len(census)))


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0], int(a[1]), int(a[2]), int(a[3]), int(a[4]), a[5])
