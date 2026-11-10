"""
The configuration file, the GNU parallel runner and the pipeline's markers.
"""
import os
import shutil
import filecmp

import pytest

from k3fq.basic import basic_stage
from k3fq.config import Config, DEFAULT_CONFIG, load_config
from k3fq.parallel import Runner
from k3fq.pipeline import Pipeline, STAGES

DATA = os.path.join(os.path.dirname(__file__), "data", "q2_sample.txt")


def write_config(tmp_path, extra=""):
    fname = tmp_path / "test.toml"
    fname.write_text('data = "%s"\njobs = 2\nchunk_size = 7\n%s\n'
                     '[[field]]\nq = 2\nmin_geometric_rank = 0\n\n'
                     '[[field]]\nq = 4\nmin_geometric_rank = 18\nclasses = 3\n' % (tmp_path / "data", extra))
    return load_config(str(fname))


def test_default_config():
    cfg = load_config(DEFAULT_CONFIG)
    assert cfg.datadir == "/scratch/k3"
    assert cfg.polred == "polredabs"
    assert cfg.qs == [2, 3, 4, 5, 7, 8]
    assert sum(f.classes for f in cfg.fields) == 13346626
    assert cfg.field(3).polyfile_name == "q3_rho16.txt"


def test_config_checks():
    with pytest.raises(ValueError):
        Config({"data": "x", "field": [{"q": 6}]})
    with pytest.raises(ValueError):
        Config({"data": "x", "field": [{"q": 2, "min_geometric_rank": 3}]})
    with pytest.raises(ValueError):
        # the cutoffs must not decrease from F_2 to F_4
        Config({"data": "x", "field": [{"q": 2, "min_geometric_rank": 20}, {"q": 4, "min_geometric_rank": 18}]})
    with pytest.raises(ValueError):
        Config({"data": "x", "number_fields": {"polred": "polredfast"}})


def sample_polyfile(tmp_path, n=20):
    polyfile = tmp_path / "q2_sample.txt"
    with open(DATA) as F, open(polyfile, "w") as G:
        for _ in range(n):
            G.write(F.readline())
    return str(polyfile)


@pytest.mark.skipif(shutil.which("parallel") is None or shutil.which("sage") is None,
                    reason="needs GNU parallel and sage")
def test_gnu_parallel_matches_in_process(tmp_path):
    cfg = write_config(tmp_path)
    polyfile = sample_polyfile(tmp_path)
    assert basic_stage(polyfile, cfg.datadir, Runner(cfg), chunksize=7) == 20
    local = str(tmp_path / "local")
    basic_stage(polyfile, local, Runner(), chunksize=7)
    names = sorted(os.listdir(os.path.join(local, "basic", "q2")))
    assert names == ["q2_sample.%05d.txt" % i for i in range(3)]
    assert filecmp.cmpfiles(os.path.join(local, "basic", "q2"), os.path.join(cfg.datadir, "basic", "q2"),
                            names, shallow=False)[0] == names
    # the joblog records three successful jobs
    logs = [f for f in os.listdir(os.path.join(cfg.datadir, "jobs")) if f.endswith(".log")]
    with open(os.path.join(cfg.datadir, "jobs", logs[0])) as F:
        assert [line.split("\t")[6] for line in F][1:] == ["0"] * 3


def test_basic_stage_rejects_stale_chunks(tmp_path):
    polyfile = sample_polyfile(tmp_path)
    datadir = str(tmp_path / "data")
    basic_stage(polyfile, datadir, chunksize=7)
    with pytest.raises(RuntimeError):
        basic_stage(polyfile, datadir, chunksize=8)


def test_rerun_markers(tmp_path):
    cfg = write_config(tmp_path)
    pipeline = Pipeline(cfg)
    for stage in STAGES:
        for st, q in pipeline._steps(stage):
            pipeline._mark_done(st, q, "")
    assert [q for st, q in pipeline._steps("census")] == [2]
    pipeline.rerun("fields", [4])
    assert not pipeline.is_done("fields", 4) and pipeline.is_done("fields", 2)
    assert not pipeline.is_done("export", 4) and pipeline.is_done("export", 2)
    assert not pipeline.is_done("combine", None)
    assert pipeline.is_done("basechange", None) and pipeline.is_done("twists", 4)
    pipeline.rerun("basic")
    assert not any(pipeline.is_done(st, q) for stage in STAGES[1:] for st, q in pipeline._steps(stage))
    assert pipeline.is_done("enumerate", 2)


def test_count_check(tmp_path):
    cfg = write_config(tmp_path)
    pipeline = Pipeline(cfg)
    assert pipeline._check("x", 2, 17) == "17 classes"
    assert pipeline._check("x", 4, 3) == "3 classes"
    with pytest.raises(RuntimeError):
        pipeline._check("x", 4, 4)
