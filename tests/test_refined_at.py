import os
import random

from sage.all import ZZ, PolynomialRing

from k3fq.basic import basic_stage
from k3fq.invariants import K3Poly
from k3fq.labels import label_to_poly
from k3fq.pgformat import read_file
from k3fq.refined_at import DELTA, base_change_orders, refined_artin_tate_failures
from k3fq.stages import QUARTIC_CENSUS_F2, export_stage, quartic_census_labels, refined_at_stage

R = PolynomialRing(ZZ, "t")
t = R.gen()

# classes over F_2 (geometric Picard rank 4) that pass the basic conditions
FAILING = [
    [2, 1, 0, 0, -2, -1, 1, 1, 1, 1, 0, 0, 0, -1, -1, -1, -1, 1, 2, 0, 0, -1, -2],
    [2, 5, 2, -4, -4, -1, 1, 3, 3, 1, 0, 0, 0, -1, -3, -3, -1, 1, 4, 4, -2, -5, -2],
    [2, 1, 0, 2, 1, 0, 1, 0, 0, 1, 0, 0, 0, -1, 0, 0, -1, 0, -1, -2, 0, -1, -2],
    [2, 3, -1, -3, -1, 0, 1, 1, -1, -1, 0, 0, 0, -1, -1, 1, 1, 0, -1, -3, -1, 3, 2],
]


def failures(coeffs, first=True):
    orders, C, S = K3Poly(coeffs)._split()
    return refined_artin_tate_failures(orders, S.list(), first=first)


def test_base_change_orders():
    assert base_change_orders([1, 2, 4, 3], 2) == [1, 1, 2, 2, 3]
    assert base_change_orders([6], 3) == [2, 2]
    assert base_change_orders([5], 5) == [1] * 4


def test_delta_table():
    assert DELTA[1] == 4 and DELTA[8] == 256 and DELTA[9] == 512 and DELTA[16] == 256
    assert set(DELTA) == set(range(1, 24))


def test_failing_examples():
    for coeffs in FAILING:
        assert K3Poly(coeffs).compute()["realizable_k3"] == 0
        assert failures(coeffs)
    assert failures(FAILING[0], first=False) == [(1, 4), (2, 2)]


def test_quadratic_case():
    # rho = 1 and three factors 1 + t: Auel-Shaikh Theorem 1 says 2 P_2(-1) >= delta_3 = 32
    P = R(FAILING[1])
    P2 = P // ((1 - t) * (1 + t) ** 3)
    assert P2 * (1 - t) * (1 + t) ** 3 == P
    assert 2 * P2(-1) < 32
    assert failures(FAILING[1], first=False) == [(1, 2)]


def test_census_sample_passes():
    labels = sorted(quartic_census_labels())
    for lab in random.Random(1).sample(labels, 300):
        assert not failures(label_to_poly(lab)), lab


def test_stage_and_export(tmp_path):
    polyfile = tmp_path / "q2_test.txt"
    polys = FAILING + [label_to_poly(lab) for lab in sorted(quartic_census_labels())[:3]]
    with open(polyfile, "w") as F:
        for c in polys:
            F.write(" ".join(map(str, c)) + "\n")
    datadir = str(tmp_path / "data")
    basic_stage(str(polyfile), datadir)
    assert refined_at_stage(2, datadir) == 4
    rows = {r["label"]: r for r in read_file(export_stage(2, datadir))}
    assert len(rows) == 7
    assert sorted(r["realizable_k3"] for r in rows.values()) == [-1] * 4 + [0] * 3
    assert all(r["realizable_kummer"] is None for r in rows.values())
