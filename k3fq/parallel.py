"""
Running worker tasks, with GNU parallel or in this process.

A worker is a function ``f(cfg, *args)`` registered with ``@worker(name)``,
whose arguments are strings (file names and numbers) and whose results are
written to files.  ``Runner.run(name, tasks)`` runs it once for each tuple of
arguments in ``tasks``:

- with GNU parallel, each task in a fresh ``sage -python k3fq.py --config
  FILE worker NAME ARGS...`` process, stopping at the first failure;
- or one after another in this process, when there is no configuration file
  (as in the tests) or only one job.
"""
import os
import subprocess
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "k3fq.py")

WORKERS = {}


def worker(name):
    def register(f):
        WORKERS[name] = f
        return f
    return register


def load_workers():
    """
    Import the modules that define workers.
    """
    from . import basic, stages  # noqa: F401


class Runner:
    def __init__(self, cfg=None):
        self.cfg = cfg

    def run(self, name, tasks, jobs=None):
        tasks = [tuple(str(a) for a in t) for t in tasks]
        if not tasks:
            return
        load_workers()
        cfg = self.cfg
        if jobs is None:
            jobs = cfg.jobs if cfg is not None else 1
        jobs = max(1, min(jobs, len(tasks)))
        if cfg is None or cfg.path is None or jobs == 1:
            for t in tasks:
                WORKERS[name](cfg, *t)
            return
        ncols = len(tasks[0])
        if any(len(t) != ncols for t in tasks) or any("\t" in a or "\n" in a for t in tasks for a in t):
            raise ValueError("bad task arguments for %s" % name)
        jobdir = os.path.join(cfg.datadir, "jobs")
        os.makedirs(jobdir, exist_ok=True)
        stamp = time.strftime("%Y%m%d-%H%M%S")
        argfile = os.path.join(jobdir, "%s.%s.args" % (name, stamp))
        with open(argfile, "w") as F:
            for t in tasks:
                F.write("\t".join(t) + "\n")
        # --quote: parallel runs the command through a shell; keep paths with spaces whole
        cmd = [cfg.parallel, "--will-cite", "--quote", "--jobs", str(jobs), "--halt", "now,fail=1",
               "--joblog", os.path.join(jobdir, "%s.%s.log" % (name, stamp)),
               "--colsep", "\t", "-a", argfile,
               cfg.sage, "-python", SCRIPT, "--config", cfg.path, "worker", name]
        cmd += ["{%d}" % (i + 1) for i in range(ncols)]
        t0 = time.time()
        print("%s: %s tasks with GNU parallel (%s jobs)" % (name, len(tasks), jobs), flush=True)
        res = subprocess.run(cmd)
        if res.returncode:
            raise RuntimeError("GNU parallel failed for %s (exit %s); see %s"
                               % (name, res.returncode, cmd[cmd.index("--joblog") + 1]))
        print("%s: done in %.1fs" % (name, time.time() - t0), flush=True)
