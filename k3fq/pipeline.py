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

A finished stage leaves a marker in done/, and later runs skip it, so an
interrupted run resumes where it stopped (the basic, fields and twists stages
also reuse their finished chunks).  ``rerun`` removes the markers of a stage
and of the stages that use its output.
"""
import os
import time

from .parallel import Runner
from .pgformat import count_rows, read_header

STAGES = ["enumerate", "basic", "fields", "twists", "kummer", "census", "basechange", "export", "combine"]
# the stages that use the output of each stage, directly or not
DEPENDENTS = {"enumerate": STAGES[1:], "basic": STAGES[2:]}
for _stage in STAGES[2:]:
    DEPENDENTS[_stage] = ["export", "combine"] if _stage not in ("export", "combine") else STAGES[STAGES.index(_stage) + 1:]


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

    def is_done(self, stage, q):
        return os.path.exists(self._marker(stage, q))

    def _mark_done(self, stage, q, note):
        os.makedirs(self.donedir, exist_ok=True)
        with open(self._marker(stage, q), "w") as F:
            F.write("%s\t%s\n" % (time.strftime("%Y-%m-%d %H:%M:%S"), note))

    def status(self):
        """
        Print a table of the finished stages (the markers in done/ also
        record when and how long they took).
        """
        qs = self.cfg.qs
        print("data directory %s" % self.datadir)
        print(("%-12s" % "" + "".join("%-8s" % ("q=%s" % q) for q in qs)).rstrip())
        for stage in STAGES:
            steps = self._steps(stage)
            if steps[0][1] is None:
                cells = ["done" if self.is_done(stage, None) else "-"] + [""] * (len(qs) - 1)
            else:
                done = {q: self.is_done(stage, q) for st, q in steps}
                cells = [("done" if done[q] else "-") if q in done else "" for q in qs]
            print("%-12s" % stage + "".join("%-8s" % c for c in cells).rstrip())

    def rerun(self, stage, qs=None):
        """
        Remove the markers of ``stage`` (for the fields ``qs``, default all)
        and of the stages that use its output.
        """
        if stage not in STAGES:
            raise ValueError("unknown stage %s; the stages are %s" % (stage, ", ".join(STAGES)))
        qs = set(self.cfg.qs if not qs else qs)
        for st in [stage] + DEPENDENTS[stage]:
            for s, q in self._steps(st):
                if q is None or q in qs:
                    marker = self._marker(s, q)
                    if os.path.exists(marker):
                        os.remove(marker)
                        print("removed %s" % marker)

    # ------------------------------------------------------------------

    def run(self, stages=None):
        os.makedirs(self.datadir, exist_ok=True)
        print("data directory %s; fields %s" % (self.datadir, ", ".join(
            "q=%s (rank >= %s)" % (f.q, f.min_geometric_rank) for f in self.cfg.fields)), flush=True)
        for stage in stages or STAGES:
            for st, q in self._steps(stage):
                if self.is_done(st, q):
                    continue
                what = st if q is None else "%s q=%s" % (st, q)
                print("\n=== %s" % what, flush=True)
                t0 = time.time()
                note = getattr(self, "_" + st)(q) or ""
                elapsed = time.time() - t0
                self._mark_done(st, q, "%.0fs %s" % (elapsed, note))
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
        self.runner.run("export", [(q0, self.datadir) for q0 in todo], jobs=self.cfg.export_jobs)
        for q0 in todo:
            n = count_rows(self._export_file(q0))
            note = self._check("export", q0, n)
            if q0 != q:
                self._mark_done("export", q0, note)
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
