import os

from k3fq.basic import basic_stage
from k3fq.invariants import K3Poly
from k3fq.labels import label_to_poly
from k3fq.pgformat import read_file
from k3fq.stages import QUARTIC_CENSUS_F2, quartic_census_labels, census_stage, export_stage

DATA = os.path.join(os.path.dirname(__file__), "data", "q2_sample.txt")


def census_head(tmp_path, n):
    fname = tmp_path / "census.txt"
    with open(QUARTIC_CENSUS_F2) as F, open(fname, "w") as G:
        for _ in range(n):
            G.write(F.readline())
    return str(fname)


def test_census_labels(tmp_path):
    labels = quartic_census_labels(census_head(tmp_path, 40))
    assert len(labels) == 40
    for lab in labels:
        D = K3Poly(label_to_poly(lab)).compute()
        # a smooth quartic has Picard rank >= 1 and passes every necessary condition
        assert D["picard_rank"] >= 1
        assert D["realizable_k3"] == 0


def test_census_stage(tmp_path):
    census = census_head(tmp_path, 5)
    polys = [label_to_poly(lab) for lab in sorted(quartic_census_labels(census))]
    with open(DATA) as F:
        others = [[int(c) for c in line.split()] for line in F][:5]
    polyfile = tmp_path / "q2_test.txt"
    with open(polyfile, "w") as F:
        for c in polys + others:
            F.write(" ".join(map(str, c)) + "\n")
    datadir = str(tmp_path / "data")
    basic_stage(str(polyfile), datadir)
    assert census_stage(2, datadir, census_file=census) == 5
    rows = {r["label"]: r for r in read_file(export_stage(2, datadir))}
    known = quartic_census_labels(census)
    for lab, r in rows.items():
        if lab in known:
            assert r["realizable_k3"] == 1
            assert r["polarization_degrees"] == [4] and r["polarization_degrees_complete"] is False
        else:
            assert r["realizable_k3"] in (0, -1)
            assert r["polarization_degrees"] is None and r["polarization_degrees_complete"] is None
    assert sum(1 for lab in rows if lab in known) == 5
