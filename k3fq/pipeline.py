"""
The pipeline that builds the k3_fq_isog upload file from a configuration
(see config.py).  Stages, in order:

    enumerate   polys/q{q}_rho{R}.txt: the polynomials with geometric Picard
                rank at least R (checked against the expected count)
    basic       basic/q{q}/: invariants of each polynomial
    fields      fields/q{q}.txt: Galois groups (Magma), reduced polynomials
                (polredabs) and LMFDB labels of the transcendental factors
    twists      twists/q{q}.txt: twist classes and minimal twist degrees
    kummer      kummer/q{q}.txt: Kummer candidates (LMFDB abelian surfaces)
    census      census/q{q}.txt: the quartic census over F_2
    basechange  basechange/q{q}.txt: primitive models (all fields at once)
    export      export/k3_fq_isog.q{q}.txt
    combine     k3_fq_isog.txt: the upload file for all fields

A finished step (a stage, for one field or for all of them) leaves a marker
in done/, recording when it finished and the parts of the configuration it
depends on.  A later run skips a step whose marker is current: present, made
with the same configuration, and newer than the markers of the steps whose
output it reads (as in make).  So an interrupted run resumes where it stopped,
and after a change of configuration, or a ``rerun`` (which removes markers),
exactly the affected steps are redone.  The basic, fields and twists stages
also reuse their finished chunks.

Markers written before 1649994 are plain text and do not record their
configuration.  ``run`` stops when it finds one; ``adopt_markers`` converts
them, for when the configuration has not changed since they were written.
"""
import json
import os
import time

from .parallel import Runner
from .pgformat import count_rows, read_header

STAGES = ["enumerate", "basic", "fields", "twists", "kummer", "census", "basechange", "export", "combine"]


class Pipeline:
    def __init__(self, cfg):
        self.cfg = cfg
        self.datadir = cfg.datadir
        self.runner = Runner(cfg)
        self.donedir = os.path.join(self.datadir, "done")

    # ------------------------------------------------------------------
    # markers

    def _steps(self, stage):
        if stage in ("basechange", "combine"):
            return [(stage, None)]
        if stage == "census":
            # only the quartic census over F_2 so far
            return [(stage, q) for q in self.cfg.qs if q == 2]
        return [(stage, q) for q in self.cfg.qs]

    def _marker(self, stage, q):
        return os.path.join(self.donedir, stage if q is None else "%s.q%s" % (stage, q))

    def _inputs(self, stage, q):
        """
        The steps whose output the step (stage, q) reads.
        """
        qs = self.cfg.qs
        if stage == "enumerate":
            return []
        if stage == "basic":
            return [("enumerate", q)]
        if stage in ("fields", "twists", "kummer", "census"):
            return [("basic", q)]
        if stage == "basechange":
            return [("basic", q0) for q0 in qs]
        if stage == "export":
            return [step for st in ("basic", "fields", "twists", "kummer", "census")
                    for step in self._steps(st) if step[1] == q] + [("basechange", None)]
        if stage == "combine":
            return [("export", q0) for q0 in qs]
        raise ValueError(stage)

    def _key(self, stage, q):
        """
        The parts of the configuration that the output of (stage, q) depends
        on, apart from its inputs.
        """
        cfg = self.cfg
        if q is None:
            # basechange and combine: the whole list of fields
            key = {"fields": [[f.q, f.min_geometric_rank, f.classes] for f in cfg.fields]}
        else:
            f = cfg.field(q)
            key = {"field": [f.q, f.min_geometric_rank, f.classes]}
        if stage == "fields":
            key.update(polred=cfg.polred,
                       polredabs_timeout=cfg.polredabs_timeout if cfg.polred == "polredabs" else None,
                       lmfdb=cfg.lmfdb_path if cfg.lmfdb_labels else None)
        elif stage == "kummer":
            key.update(lmfdb=cfg.lmfdb_path)
        elif stage == "census":
            key.update(census=cfg.census_file)
        return json.loads(json.dumps(key))

    def _read_marker(self, stage, q):
        fname = self._marker(stage, q)
        if not os.path.exists(fname):
            return None
        with open(fname) as F:
            text = F.read()
        try:
            return json.loads(text)
        except ValueError:
            # written before markers recorded their configuration
            return {"legacy": True, "time": os.path.getmtime(fname), "note": text.strip()}

    def state(self, stage, q):
        """
        "done", "-" (never finished), "config" (finished with another
        configuration), "inputs" (an input has changed since), or "legacy"
        (an old marker, made with an unknown configuration).
        """
        m = self._read_marker(stage, q)
        if m is None:
            return "-"
        if m.get("legacy"):
            return "legacy"
        if m["key"] != self._key(stage, q):
            return "config"
        for st, q0 in self._inputs(stage, q):
            if self.state(st, q0) != "done" or self._read_marker(st, q0)["time"] > m["time"]:
                return "inputs"
        return "done"

    def is_done(self, stage, q):
        return self.state(stage, q) == "done"

    def _mark_done(self, stage, q, seconds, note, when=None):
        os.makedirs(self.donedir, exist_ok=True)
        when = time.time() if when is None else when
        marker = {"time": when, "finished": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(when)),
                  "seconds": round(seconds), "note": note, "key": self._key(stage, q)}
        fname = self._marker(stage, q)
        with open(fname + ".tmp", "w") as F:
            json.dump(marker, F)
            F.write("\n")
        os.replace(fname + ".tmp", fname)

    def status(self):
        """
        Print a table of the steps: done, - (not yet), or redo (the
        configuration or an input changed since it finished).
        """
        qs = self.cfg.qs
        print("data directory %s" % self.datadir)
        print(("%-12s" % "" + "".join("%-8s" % ("q=%s" % q) for q in qs)).rstrip())
        show = {"done": "done", "-": "-", "config": "redo", "inputs": "redo", "legacy": "old"}
        for stage in STAGES:
            steps = self._steps(stage)
            if steps and steps[0][1] is None:
                cells = [show[self.state(stage, None)]]
            else:
                states = {q: show[self.state(stage, q)] for st, q in steps}
                cells = [states.get(q, "") for q in qs]
            print(("%-12s" % stage + "".join("%-8s" % c for c in cells)).rstrip())

    def legacy_markers(self):
        return [(st, q) for stage in STAGES for st, q in self._steps(stage)
                if (self._read_marker(st, q) or {}).get("legacy")]

    def adopt_markers(self):
        """
        Convert the old plain-text markers of the configured steps, recording
        the current configuration in them (keeping their times).  Only right
        if the configuration has not changed since they were written.
        """
        steps = self.legacy_markers()
        for st, q in steps:
            m = self._read_marker(st, q)
            self._mark_done(st, q, 0, "adopted: %s" % m["note"], when=m["time"])
            print("adopted %s" % self._marker(st, q))
        if not steps:
            print("no old markers")

    def rerun(self, stage, qs=None):
        """
        Remove the markers of ``stage`` (for the fields ``qs``, default all);
        the steps that read its output are then redone too.
        """
        if stage not in STAGES:
            raise ValueError("unknown stage %s; the stages are %s" % (stage, ", ".join(STAGES)))
        qs = set(self.cfg.qs if not qs else qs)
        for st, q in self._steps(stage):
            if q is None or q in qs:
                marker = self._marker(st, q)
                if os.path.exists(marker):
                    os.remove(marker)
                    print("removed %s" % marker)

    # ------------------------------------------------------------------

    def run(self, stages=None):
        legacy = self.legacy_markers()
        if legacy:
            raise RuntimeError(
                "%s has %s markers from before 1649994 (e.g. %s), which do not record the configuration "
                "they were made with.  If the configuration has not changed since, keep that work with "
                "`k3fq.py adopt-markers`; otherwise remove them (or the whole done/ directory) to redo it."
                % (self.donedir, len(legacy), os.path.basename(self._marker(*legacy[0]))))
        os.makedirs(self.datadir, exist_ok=True)
        print("data directory %s; fields %s" % (self.datadir, ", ".join(
            "q=%s (rank >= %s)" % (f.q, f.min_geometric_rank) for f in self.cfg.fields)), flush=True)
        reasons = {"-": "", "config": " (the configuration has changed)", "inputs": " (its input has changed)"}
        for stage in stages or STAGES:
            for st, q in self._steps(stage):
                state = self.state(st, q)
                if state == "done":
                    continue
                what = st if q is None else "%s q=%s" % (st, q)
                print("\n=== %s%s" % (what, reasons[state]), flush=True)
                t0 = time.time()
                note = getattr(self, "_" + st)(q) or ""
                elapsed = time.time() - t0
                self._mark_done(st, q, elapsed, note)
                print("=== %s done (%.0fs) %s" % (what, elapsed, note), flush=True)

    def _check(self, what, q, n):
        expected = self.cfg.field(q).classes
        if expected is not None and n != expected:
            raise RuntimeError("%s for q=%s: %s classes, but %s were expected" % (what, q, n, expected))
        return "%s classes" % n

    def _polyfile(self, q):
        return os.path.join(self.datadir, "polys", self.cfg.field(q).polyfile_name)

    # ------------------------------------------------------------------
    # stages

    def _enumerate(self, q):
        from .enumerate import enumerate_elevated
        f = self.cfg.field(q)
        n = enumerate_elevated(q, f.min_geometric_rank, os.path.dirname(self._polyfile(q)),
                               threads=self.cfg.jobs)
        return self._check("enumeration", q, n)

    def _basic(self, q):
        from .basic import basic_stage
        from .stages import basic_files
        basic_stage(self._polyfile(q), self.datadir, self.runner, self.cfg.chunk_size)
        n = sum(count_rows(f) for f in basic_files(self.datadir, q))
        return self._check("basic stage", q, n)

    def _fields(self, q):
        from .stages import fields_stage
        return "%s factors" % fields_stage(q, self.datadir, self.runner)

    def _twists(self, q):
        from .stages import twists_stage
        return self._check("twists stage", q, twists_stage(q, self.datadir, self.runner))

    def _kummer(self, q):
        from .stages import kummer_stage
        kummer_stage(q, self.datadir, self.cfg)

    def _census(self, q):
        from .stages import census_stage
        return "%s classes" % census_stage(q, self.datadir, self.cfg.census_file)

    def _basechange(self, q):
        from .stages import basechange_stage
        basechange_stage(self.datadir, self.cfg.qs)

    def _export(self, q):
        # all fields at once, in parallel; the markers record each field
        todo = [q0 for q0 in self.cfg.qs if not self.is_done("export", q0)]
        t0 = time.time()
        self.runner.run("export", [(q0, self.datadir) for q0 in todo], jobs=self.cfg.export_jobs)
        for q0 in todo:
            note = self._check("export", q0, count_rows(self._export_file(q0)))
            if q0 != q:
                self._mark_done("export", q0, time.time() - t0, note)
        return self._check("export", q, count_rows(self._export_file(q)))

    def _export_file(self, q):
        return os.path.join(self.datadir, "export", "k3_fq_isog.q%s.txt" % q)

    def _combine(self, q):
        files = [self._export_file(q0) for q0 in self.cfg.qs]
        header = read_header(files[0])
        for f in files[1:]:
            if read_header(f) != header:
                raise RuntimeError("%s has different columns from %s" % (f, files[0]))
        out = os.path.join(self.datadir, "k3_fq_isog.txt")
        tmp = out + ".tmp"
        total = 0
        with open(tmp, "wb") as F:
            for i, fname in enumerate(files):
                with open(fname, "rb") as G:
                    for _ in range(3):
                        line = G.readline()
                        if i == 0:
                            F.write(line)
                    while True:
                        block = G.read(1 << 24)
                        if not block:
                            break
                        total += block.count(b"\n")
                        F.write(block)
        expected = [f.classes for f in self.cfg.fields]
        if None not in expected and total != sum(expected):
            raise RuntimeError("%s has %s rows, but %s were expected" % (out, total, sum(expected)))
        os.replace(tmp, out)
        print("wrote %s" % out)
        return "%s classes" % total
