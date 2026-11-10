"""
Regression tests for issues found in review (2026-11-05).
"""
import os

from k3fq.cyclotomic import R, T
from k3fq.invariants import K3Poly
from k3fq.columns import BASIC_COLUMNS, ISOG_TYPES
from k3fq.pgformat import write_file, read_file


def test_ordinary_needs_slope_zero():
    # q = 4, P = 2(2 - T + 2T^2)(1 - T)^20: slopes 1/2, 1 (x20), 3/2
    D = K3Poly((2 * R([2, -1, 2]) * (1 - T) ** 20).list()).compute()
    assert D["poly"][1] == -82
    assert D["height"] == 1 and not D["has_k3_newton_polygon"]
    assert not D["is_ordinary"]


def test_wreath_table_used():
    # degree 12 factor whose full-wreath certificate needs data/wreath_labels.txt
    f = R([2, -1, 2, -1, 3, -1, 2, -1, 3, -1, 2, -1, 2])
    D = K3Poly((f * (1 - T) ** 10).list()).compute()
    assert D["galois_groups"][0] is not None
    assert D["galois_groups"][0].startswith("12T")


def test_twist_stage_deduplicates(tmp_path):
    from k3fq.stages import twists_stage
    path = os.path.join(os.path.dirname(__file__), "data", "q2_sample.txt")
    with open(path) as F:
        coeffs = [int(c) for c in F.readline().split()]
    row = K3Poly(coeffs).compute()
    d = tmp_path / "basic" / "q2"
    d.mkdir(parents=True)
    cols = [(c, ISOG_TYPES[c]) for c in BASIC_COLUMNS]
    write_file(str(d / "a.00000.txt"), cols, [row])
    write_file(str(d / "b.00000.txt"), cols, [row])
    twists_stage(2, str(tmp_path))
    out = list(read_file(str(tmp_path / "twists" / "q2.txt")))
    assert len(out) == 1
    assert out[0]["twist_count"] == 1 and out[0]["twists"] == []
