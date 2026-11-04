import os

from sage.all import ZZ

from k3fq.cyclotomic import R, T
from k3fq.invariants import K3Poly
from k3fq.twists import base_change, base_change_resultant, twist_data, primitive_model_data

DATA = os.path.join(os.path.dirname(__file__), "data", "q2_sample.txt")


def sample(n):
    with open(DATA) as F:
        return [[int(c) for c in line.split()] for line in F][:n]


def test_base_change_composition():
    for c in sample(5):
        P = R(c)
        P2 = base_change(c, 2)
        PP = P * P(-T)
        assert R(PP.list()[::2]) == R(P2)
        assert base_change(c, 6) == base_change(base_change(c, 2), 3)
        assert R(base_change(c, 3))[0] == 8
        for n in (2, 3, 5, 6):
            assert base_change(c, n) == base_change_resultant(c, n)


def test_quadratic_twists():
    rows = [K3Poly(c).compute() for c in sample(15)]
    twisted = [K3Poly(R(r["poly"])(-T).list()).compute() for r in rows]
    td = twist_data(rows + twisted)
    for r, s in zip(rows, twisted):
        if r["label"] == s["label"]:
            continue
        a, b = td[r["label"]], td[s["label"]]
        assert a["geometric_class"] == b["geometric_class"]
        degs = [d for lab, bc, d in a["twists"] if lab == s["label"]]
        assert 2 in degs or 1 in degs


def test_primitive_models():
    rows2 = [K3Poly(c).compute() for c in sample(10)]
    rows4 = [K3Poly(base_change(r["poly"], 2)).compute() for r in rows2]
    data = primitive_model_data({2: rows2, 4: rows4})
    for r2, r4 in zip(rows2, rows4):
        assert data[r2["label"]]["is_primitive"]
        assert r2["label"] in data[r4["label"]]["primitive_models"]
        assert data[r4["label"]]["is_primitive"] is False
