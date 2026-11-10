"""
The configuration file of the pipeline, in TOML (see k3fq.toml).

    data = "/scratch/k3"          # data directory (all intermediate files and the output)
    jobs = 0                      # simultaneous jobs; 0 means one per core
    chunk_size = 20000            # classes per job in the basic stage
    export_jobs = 1               # simultaneous exports (each holds a field's data in memory)

    [programs]
    parallel = "parallel"         # GNU parallel
    sage = "sage"                 # workers run as `sage -python k3fq.py worker ...`
    magma = "magma"
    magma_jobs = 1                # simultaneous Magma processes

    [number_fields]
    polred = "polredabs"          # or "polredbest"
    polredabs_timeout = 600       # seconds per polynomial before falling back to polredbest
    lmfdb_labels = true           # look up the reduced polynomials in the LMFDB's nf_fields

    [lmfdb]
    path = "~roed/lmfdb"          # an LMFDB checkout, for av_fq_isog and nf_fields

    [census]
    quartics_f2 = "external/root-unitary/k3-quartic-f2/k3final.txt"

    [[field]]                     # one per base field
    q = 2
    min_geometric_rank = 0
    classes = 5478058             # expected number of classes (optional, checked)
"""
import os
import tomllib

from sage.all import ZZ

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_CONFIG = os.path.join(ROOT, "k3fq.toml")


class Field:
    def __init__(self, q, min_geometric_rank, classes=None):
        self.q = int(q)
        self.min_geometric_rank = int(min_geometric_rank)
        self.classes = None if classes is None else int(classes)

    @property
    def polyfile_name(self):
        return "q%s_rho%s.txt" % (self.q, self.min_geometric_rank)


class Config:
    def __init__(self, data, path=None):
        self.path = None if path is None else os.path.abspath(path)
        self.datadir = os.path.abspath(os.path.expanduser(data["data"]))
        jobs = int(data.get("jobs", 0))
        self.jobs = jobs if jobs > 0 else os.cpu_count()
        self.chunk_size = int(data.get("chunk_size", 20000))
        self.export_jobs = int(data.get("export_jobs", 1))
        programs = data.get("programs", {})
        self.parallel = programs.get("parallel", "parallel")
        self.sage = programs.get("sage", "sage")
        self.magma = programs.get("magma", "magma")
        self.magma_jobs = int(programs.get("magma_jobs", 1))
        nf = data.get("number_fields", {})
        self.polred = nf.get("polred", "polredabs")
        if self.polred not in ("polredabs", "polredbest"):
            raise ValueError("number_fields.polred must be polredabs or polredbest")
        self.polredabs_timeout = float(nf.get("polredabs_timeout", 600))
        self.lmfdb_labels = bool(nf.get("lmfdb_labels", False))
        self.lmfdb_path = data.get("lmfdb", {}).get("path")
        if self.lmfdb_path:
            self.lmfdb_path = os.path.expanduser(self.lmfdb_path)
        census = data.get("census", {}).get("quartics_f2", os.path.join(
            "external", "root-unitary", "k3-quartic-f2", "k3final.txt"))
        self.census_file = census if os.path.isabs(census) else os.path.join(ROOT, census)
        self.fields = [Field(f["q"], f.get("min_geometric_rank", 0), f.get("classes"))
                       for f in data.get("field", [])]
        self._check_fields()

    def _check_fields(self):
        by_q = {}
        for f in self.fields:
            if not ZZ(f.q).is_prime_power():
                raise ValueError("q = %s is not a prime power" % f.q)
            if f.q in by_q:
                raise ValueError("q = %s appears twice" % f.q)
            if f.min_geometric_rank % 2 or not 0 <= f.min_geometric_rank <= 22:
                raise ValueError("the minimal geometric Picard rank must be even and at most 22")
            by_q[f.q] = f
        # The primitivity logic needs the cutoffs to be nondecreasing along subfields
        for f in self.fields:
            p, r = ZZ(f.q).is_prime_power(get_data=True)
            for d in ZZ(r).divisors():
                q0 = int(p ** (r // d))
                if q0 != f.q and q0 in by_q and by_q[q0].min_geometric_rank > f.min_geometric_rank:
                    raise ValueError("the cutoff for q = %s is larger than for q = %s" % (q0, f.q))

    @property
    def qs(self):
        return [f.q for f in self.fields]

    def field(self, q):
        for f in self.fields:
            if f.q == q:
                return f
        raise KeyError("q = %s is not in the configuration" % q)


def load_config(path=None):
    path = path or DEFAULT_CONFIG
    with open(path, "rb") as F:
        return Config(tomllib.load(F), path)
