"""
The configuration file, the GNU parallel runner and the pipeline's markers.
"""
import os
import shutil
import sys
import filecmp

import pytest

from k3fq.basic import basic_stage
from k3fq.config import Config, DEFAULT_CONFIG, load_config
from k3fq.parallel import Runner
from k3fq.pipeline import Pipeline, STAGES

DATA = os.path.join(os.path.dirname(__file__), "data", "q2_sample.txt")


FIELDS = ('[[field]]\nq = 2\nmin_geometric_rank = 0\n\n'
          '[[field]]\nq = 4\nmin_geometric_rank = 18\nclasses = 3\n')


def write_config(tmp_path, extra="", fields=FIELDS):
    # spaces in the paths, which GNU parallel must keep whole
    fname = tmp_path / "test config.toml"
    fname.write_text('data = "%s"\njobs = 2\nchunk_size = 7\n%s\n%s'
                     % (tmp_path / "data dir", extra, fields))
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


def mark_all(pipeline):
    for stage in STAGES:
        for st, q in pipeline._steps(stage):
            pipeline._mark_done(st, q, 0, "")


def states(pipeline):
    return {(st, q): pipeline.state(st, q) for stage in STAGES for st, q in pipeline._steps(stage)}


def test_rerun_markers(tmp_path):
    pipeline = Pipeline(write_config(tmp_path))
    mark_all(pipeline)
    assert set(states(pipeline).values()) == {"done"}
    assert [q for st, q in pipeline._steps("census")] == [2]
    pipeline.rerun("fields", [4])
    S = states(pipeline)
    assert S["fields", 4] == "-" and S["fields", 2] == "done"
    assert S["export", 4] == "inputs" and S["export", 2] == "done"
    assert S["combine", None] == "inputs"
    assert S["basechange", None] == "done" and S["twists", 4] == "done"


def test_rerun_through_basechange(tmp_path):
    # base change over F_4 uses the classes over F_2, so redoing the basic
    # stage for q = 2 redoes the export for q = 4 too
    pipeline = Pipeline(write_config(tmp_path))
    mark_all(pipeline)
    pipeline.rerun("basic", [2])
    S = states(pipeline)
    assert S["basic", 2] == "-" and S["basic", 4] == "done"
    assert S["fields", 2] == "inputs" and S["fields", 4] == "done"
    assert S["basechange", None] == "inputs"
    assert S["export", 2] == S["export", 4] == S["combine", None] == "inputs"
    # redoing them in order makes everything current again
    for stage in STAGES:
        for st, q in pipeline._steps(stage):
            if pipeline.state(st, q) != "done":
                pipeline._mark_done(st, q, 0, "")
    assert set(states(pipeline).values()) == {"done"}


def test_configuration_changes(tmp_path):
    only2 = '[[field]]\nq = 2\nmin_geometric_rank = 0\n'
    mark_all(Pipeline(write_config(tmp_path, fields=only2)))
    # adding a field redoes the steps for all fields together
    S = states(Pipeline(write_config(tmp_path)))
    assert S["enumerate", 2] == S["basic", 2] == S["fields", 2] == "done"
    assert S["enumerate", 4] == "-"
    assert S["basechange", None] == S["combine", None] == "config"
    assert S["export", 2] == "inputs"
    # changing a cutoff or an expected count redoes that field from the start
    mark_all(Pipeline(write_config(tmp_path)))
    S = states(Pipeline(write_config(tmp_path, fields=FIELDS.replace("classes = 3", "classes = 4"))))
    assert S["enumerate", 4] == S["basic", 4] == S["export", 4] == "config" and S["enumerate", 2] == "done"
    # the reduction method only affects the fields stage and what reads it
    S = states(Pipeline(write_config(tmp_path, '[number_fields]\npolred = "polredbest"\n')))
    assert S["fields", 2] == S["fields", 4] == "config"
    assert S["twists", 2] == "done" and S["export", 2] == "inputs"


def test_markers_from_9179df5(tmp_path):
    pipeline = Pipeline(write_config(tmp_path))
    os.makedirs(pipeline.donedir)
    for st, q in [("enumerate", 2), ("basic", 2)]:
        with open(pipeline._marker(st, q), "w") as F:
            F.write("2026-11-10 01:00:00\t77s 182408 classes\n")
    S = states(pipeline)
    assert S["enumerate", 2] == S["basic", 2] == "done" and S["fields", 2] == "-"


def test_status_without_q2(tmp_path, capsys):
    pipeline = Pipeline(write_config(tmp_path, fields='[[field]]\nq = 3\nmin_geometric_rank = 22\n'))
    pipeline.status()
    out = capsys.readouterr().out
    assert "census" in out and "basechange  -" in out


def test_polred_chunks_depend_on_settings(tmp_path):
    from k3fq.stages import _factor_chunks
    factors = [(1, 1, 1), (2, 1, 1)]
    a = _factor_chunks(str(tmp_path), "polred", factors, "polredabs 600")
    b = _factor_chunks(str(tmp_path), "polred", factors, "polredbest")
    c = _factor_chunks(str(tmp_path), "polred", factors, "polredabs 600")
    assert a == c and a != b


def test_lmfdb_imported_from_checkout(tmp_path):
    # the LMFDB reads secrets.ini relative to the working directory
    from k3fq.lmfdb_access import lmfdb_db
    checkout = tmp_path / "lmfdb checkout"
    (checkout / "lmfdb").mkdir(parents=True)
    (checkout / "lmfdb" / "__init__.py").write_text("import os\ndb = os.getcwd()\n")
    cfg = write_config(tmp_path, '[lmfdb]\npath = "%s"\n' % checkout)
    cwd = os.getcwd()
    saved = sys.modules.pop("lmfdb", None)
    try:
        assert lmfdb_db(cfg) == str(checkout)
        assert os.getcwd() == cwd
    finally:
        sys.modules.pop("lmfdb", None)
        if saved is not None:
            sys.modules["lmfdb"] = saved
        sys.path.remove(str(checkout))


def test_count_check(tmp_path):
    cfg = write_config(tmp_path)
    pipeline = Pipeline(cfg)
    assert pipeline._check("x", 2, 17) == "17 classes"
    assert pipeline._check("x", 4, 3) == "3 classes"
    with pytest.raises(RuntimeError):
        pipeline._check("x", 4, 4)
