"""
The stages after ``basic``, and the final export.  All intermediate files
live under the data directory:

    polys/q{q}_*.txt          output of ``enumerate`` (one polynomial per line)
    basic/q{q}/*.txt          output of ``basic`` (BASIC_COLUMNS)
    fields/q{q}.txt           Galois groups and number fields, keyed by factor
    fields/q{q}/*             inputs and outputs of the Magma and polred jobs
    twists/q{q}/keys.*.txt    geometric keys (intermediate)
    twists/q{q}/pairs.*       inputs and outputs of the pairwise twist jobs
    twists/q{q}.txt           TWIST_COLUMNS
    basechange/q{q}.txt       BASECHANGE_COLUMNS
    kummer/q{q}.txt           KUMMER_COLUMNS
    census/q{q}.txt           CENSUS_COLUMNS, for the classes known to contain a K3 surface
    refined_at/q{q}.txt       the classes ruled out by the refined Artin-Tate conditions
    export/k3_fq_isog.q{q}.txt  all columns, ready for psycodict's copy_from
"""
import glob
import hashlib
import os
import subprocess
import time
from collections import defaultdict

from sage.all import ZZ, pari

from .columns import (ISOG_COLUMNS, ISOG_TYPES, TWIST_COLUMNS, BASECHANGE_COLUMNS,
                      KUMMER_COLUMNS, CENSUS_COLUMNS)
from .cyclotomic import R, split_cyclotomic
from .labels import label_sort_key, poly_to_label
from .lmfdb_access import lmfdb_db
from .parallel import Runner, worker
from .pgformat import read_file, write_file
from . import twists as tw

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def basic_files(datadir, q):
    files = sorted(glob.glob(os.path.join(datadir, "basic", "q%s" % q, "*.txt")))
    if not files:
        raise RuntimeError("no basic stage output for q=%s" % q)
    return files


def iter_rows(datadir, q, columns=None):
    seen = set()
    wanted = None if columns is None else list(dict.fromkeys(["label"] + list(columns)))
    for fname in basic_files(datadir, q):
        for row in read_file(fname, wanted):
            if row["label"] in seen:
                continue
            seen.add(row["label"])
            if columns is not None:
                row = {c: row[c] for c in columns}
            yield row


def available_qs(datadir):
    return sorted(int(os.path.basename(d)[1:]) for d in glob.glob(os.path.join(datadir, "basic", "q*")))


# ----------------------------------------------------------------------
# fields: Galois groups that the certificate missed (Magma), reduced
# defining polynomials (polredabs, or polredbest), and LMFDB labels

FACTOR_CHUNK = 2000


def _write_lines_atomic(fname, lines):
    tmp = fname + ".tmp"
    with open(tmp, "w") as F:
        F.writelines(line + "\n" for line in lines)
    os.replace(tmp, fname)


def _factor_chunks(workdir, kind, factors, settings=""):
    """
    Write ``factors`` in chunks to ``workdir/{kind}.{hash}.in`` and return
    the (input, output) pairs.  The names depend on the contents and on
    ``settings``, so an output from an earlier run is only reused for the same
    input and settings.
    """
    pairs = []
    for i in range(0, len(factors), FACTOR_CHUNK):
        lines = [",".join(map(str, f)) for f in factors[i:i + FACTOR_CHUNK]]
        h = hashlib.sha1("\n".join([settings] + lines).encode()).hexdigest()[:16]
        infile = os.path.join(workdir, "%s.%s.in" % (kind, h))
        if not os.path.exists(infile):
            _write_lines_atomic(infile, lines)
        pairs.append((infile, infile[:-3] + ".out"))
    return pairs


def reduced_polynomial(f, method="polredabs", timeout=600):
    """
    A reduced defining polynomial of Q[T]/(f) (coefficients, constant term
    first), and the method that produced it: polredabs unless it takes longer
    than ``timeout`` seconds (0 for no limit), and then polredbest.
    """
    g = pari(R(list(f)))
    if method == "polredabs":
        from sage.all import alarm, cancel_alarm
        from cysignals.signals import AlarmInterrupt
        try:
            if timeout:
                alarm(timeout)
            h = g.polredabs()
            return [int(c) for c in pari.Vecrev(h)], "polredabs"
        except AlarmInterrupt:
            pass
        finally:
            cancel_alarm()
    h = g.polredbest()
    return [int(c) for c in pari.Vecrev(h)], "polredbest"


def polred_settings(cfg):
    method = cfg.polred if cfg is not None else "polredabs"
    timeout = cfg.polredabs_timeout if cfg is not None else 600
    return method, timeout


@worker("polred")
def polred_chunk(cfg, infile, outfile):
    method, timeout = polred_settings(cfg)
    t0 = time.time()
    lines = []
    with open(infile) as F:
        for line in F:
            f = [int(c) for c in line.split(",")]
            g, used = reduced_polynomial(f, method, timeout)
            lines.append("%s\t%s\t%s" % (line.strip(), ",".join(map(str, g)), used))
    _write_lines_atomic(outfile, lines)
    print("%s: %s polynomials in %.1fs" % (os.path.basename(outfile), len(lines), time.time() - t0), flush=True)


@worker("magma")
def magma_chunk(cfg, infile, outfile):
    """
    Transitive labels of the Galois groups of the polynomials in ``infile``,
    computed by magma/galois_groups.m.
    """
    script = os.path.join(ROOT, "magma", "galois_groups.m")
    magma = cfg.magma if cfg is not None else "magma"
    t0 = time.time()
    tmp = outfile + ".tmp"
    subprocess.run([magma, "-b", "infile:=%s" % infile, "outfile:=%s" % tmp, script],
                   check=True, stdout=subprocess.DEVNULL)
    os.replace(tmp, outfile)
    print("%s: Magma in %.1fs" % (os.path.basename(outfile), time.time() - t0), flush=True)


def _read_pairs(pairs, parse):
    ans = {}
    for infile, outfile in pairs:
        with open(outfile) as F:
            for line in F:
                k, v = parse(line)
                ans[k] = v
    return ans


def fields_stage(q, datadir, runner=None):
    runner = runner or Runner()
    cfg = runner.cfg
    factors = set()
    missing = set()
    for row in iter_rows(datadir, q, ["transcendental_factors", "galois_groups"]):
        for (f, e), g in zip(row["transcendental_factors"], row["galois_groups"]):
            factors.add(tuple(f))
            if g is None:
                missing.add(tuple(f))
    print("q=%s: %s distinct transcendental factors, %s without a Galois group"
          % (q, len(factors), len(missing)), flush=True)
    workdir = os.path.join(datadir, "fields", "q%s" % q)
    os.makedirs(workdir, exist_ok=True)

    pairs = _factor_chunks(workdir, "magma", sorted(missing))
    runner.run("magma", [p for p in pairs if not os.path.exists(p[1])],
               jobs=cfg.magma_jobs if cfg is not None else 1)

    def parse_magma(line):
        coeffs, label = line.split()
        return tuple(int(c) for c in coeffs.split(",")), label
    galois = _read_pairs(pairs, parse_magma)
    assert set(galois) == missing

    factors = sorted(factors)
    method, timeout = polred_settings(cfg)
    pairs = _factor_chunks(workdir, "polred", factors,
                           "polredabs %s" % timeout if method == "polredabs" else method)
    runner.run("polred", [p for p in pairs if not os.path.exists(p[1])])

    def parse_polred(line):
        f, g, used = line.rstrip("\n").split("\t")
        return tuple(int(c) for c in f.split(",")), ([int(c) for c in g.split(",")], used)
    reduced = _read_pairs(pairs, parse_polred)
    assert set(reduced) == set(factors)
    print("q=%s: %s reduced with polredabs, %s with polredbest" % (
        q, sum(1 for g, u in reduced.values() if u == "polredabs"),
        sum(1 for g, u in reduced.values() if u == "polredbest")), flush=True)

    labels = {}
    if cfg is not None and cfg.lmfdb_labels:
        db = lmfdb_db(cfg)
        t0 = time.time()
        for f in factors:
            g, used = reduced[f]
            if used == "polredabs":
                lab = db.nf_fields.lucky({"coeffs": g}, "label")
                if lab:
                    labels[f] = lab
        print("q=%s: %s fields found in the LMFDB (%.1fs)" % (q, len(labels), time.time() - t0), flush=True)

    out = os.path.join(datadir, "fields", "q%s.txt" % q)
    rows = [{"factor": list(f), "galois_group": galois.get(f), "nf_label": labels.get(f),
             "nf_coeffs": reduced[f][0]} for f in factors]
    write_file(out, [("factor", "jsonb"), ("galois_group", "text"), ("nf_label", "text"),
                     ("nf_coeffs", "jsonb")], rows)
    return len(rows)


# ----------------------------------------------------------------------
# twists

_TWIST_INPUT = ["label", "q", "cyclotomic_factors", "geometric_picard_rank", "transcendental_poly",
                "transcendental_factors"]
PAIR_BATCH = 500   # twist classes per pairwise job


def _newer(output, inputs):
    return os.path.exists(output) and all(os.path.getmtime(output) >= os.path.getmtime(f) for f in inputs)


@worker("twistkeys")
def twist_keys_chunk(cfg, basicfile, outfile):
    rows = []
    for row in read_file(basicfile, ["label", "q", "geometric_picard_rank", "transcendental_factors"]):
        rows.append({"label": row["label"],
                     "key": hashlib.sha1(tw.geometric_key(row).encode()).hexdigest()})
    write_file(outfile, [("label", "text"), ("key", "text")], rows)


@worker("twistpairs")
def twist_pairs_batch(cfg, infile, outfile):
    classes = defaultdict(list)
    for row in read_file(infile):
        classes[row.pop("class")].append(row)
    rows = []
    for members in classes.values():
        cols = tw.twist_columns(members, tw.class_twists(members))
        rows.extend(dict(label=lab, **c) for lab, c in cols.items())
    write_file(outfile, [("label", "text")] + [(c, ISOG_TYPES[c]) for c in TWIST_COLUMNS], rows)


def twists_stage(q, datadir, runner=None):
    runner = runner or Runner()
    keydir = os.path.join(datadir, "twists", "q%s" % q)
    os.makedirs(keydir, exist_ok=True)
    t0 = time.time()
    tasks = []
    for fname in basic_files(datadir, q):
        outfile = os.path.join(keydir, "keys." + os.path.basename(fname))
        if not _newer(outfile, [fname]):
            tasks.append((fname, outfile))
    runner.run("twistkeys", tasks)
    print("q=%s: geometric keys computed (%.1fs)" % (q, time.time() - t0), flush=True)
    # sets: overlapping enumeration files (e.g. full and elevated rank) share labels
    classes = defaultdict(set)
    for fname in basic_files(datadir, q):
        for row in read_file(os.path.join(keydir, "keys." + os.path.basename(fname))):
            classes[row["key"]].add(row["label"])
    small = {lab for labs in classes.values() if 1 < len(labs) <= tw.MAX_TWISTS for lab in labs}
    data = {r["label"]: r for r in iter_rows(datadir, q, _TWIST_INPUT) if r["label"] in small}
    out = {}
    work = []
    for labs in classes.values():
        labs = sorted(labs, key=label_sort_key)
        if 1 < len(labs) <= tw.MAX_TWISTS:
            work.append(labs)
        else:
            out.update(tw.twist_columns([{"label": lab} for lab in labs],
                                        None if len(labs) > 1 else {labs[0]: []}))
    for old in glob.glob(os.path.join(keydir, "pairs.*")):
        os.remove(old)
    cols = [("class", "integer")] + [(c, ISOG_TYPES[c]) for c in _TWIST_INPUT]
    tasks = []
    for b in range(0, len(work), PAIR_BATCH):
        infile = os.path.join(keydir, "pairs.%05d.in" % (b // PAIR_BATCH))
        rows = [dict(data[lab], **{"class": b + i}) for i, labs in enumerate(work[b:b + PAIR_BATCH]) for lab in labs]
        write_file(infile, cols, rows)
        tasks.append((infile, infile[:-3] + ".out"))
    runner.run("twistpairs", tasks)
    for infile, outfile in tasks:
        out.update({row.pop("label"): row for row in read_file(outfile)})
    rows = [dict(label=lab, **c) for lab, c in sorted(out.items(), key=lambda kv: label_sort_key(kv[0]))]
    write_file(os.path.join(datadir, "twists", "q%s.txt" % q),
               [("label", "text")] + [(c, ISOG_TYPES[c]) for c in TWIST_COLUMNS], rows)
    print("q=%s: %s twist classes, %s with pairwise data (%.1fs)" % (q, len(classes), len(work), time.time() - t0), flush=True)
    return len(rows)


# ----------------------------------------------------------------------
# primitive models (needs the basic stage for all relevant q)

def basechange_stage(datadir, qs=None):
    """
    Primitive models.  Base change preserves the geometric Picard rank, so a
    class over a subfield F_{q0} can only be a model of a class in the data
    over F_q if its rank is at least the smallest rank over F_q; only those
    are base changed.  (Since the cutoffs are nondecreasing in q, the data
    over F_{q0} contains all of them.)  The classes over F_q are streamed.
    ``qs`` defaults to every q in the data directory.
    """
    qs = sorted(qs) if qs else available_qs(datadir)
    minrho = {q: min(r["geometric_picard_rank"] for r in iter_rows(datadir, q, ["geometric_picard_rank"]))
              for q in qs}
    os.makedirs(os.path.join(datadir, "basechange"), exist_ok=True)
    for q in qs:
        subs = tw.subfields(q)
        complete = all(q0 in qs for q0 in subs)
        sources = {q0: [r for r in iter_rows(datadir, q0, ["label", "poly", "geometric_picard_rank"])
                        if r["geometric_picard_rank"] >= minrho[q]]
                   for q0 in subs if q0 in qs}
        models = tw.primitive_models(q, sources)
        rows = (dict(label=r["label"], **tw.primitivity(models.get(r["label"]), complete))
                for r in iter_rows(datadir, q, ["label"]))
        write_file(os.path.join(datadir, "basechange", "q%s.txt" % q),
                   [("label", "text")] + [(c, ISOG_TYPES[c]) for c in BASECHANGE_COLUMNS], rows)


# ----------------------------------------------------------------------
# Kummer surfaces

def wedge_square_transcendental(q, a1, a2):
    """
    For an abelian surface with L-polynomial 1 + a1 T + a2 T^2 + q a1 T^3 +
    q^2 T^4, the normalized wedge square is (1-T)^2 W(T) with
    q W(T) = q - (a2 - 2q) T + (2q + a1^2 - 2 a2) T^2 - (a2 - 2q) T^3 + q T^4.
    Returns the coefficient list of the transcendental part of q W.
    """
    W = R([q, -(a2 - 2 * q), 2 * q + a1 ** 2 - 2 * a2, -(a2 - 2 * q), q])
    return split_cyclotomic(W)[2].list()


def lmfdb_abelian_surfaces(q, cfg):
    """
    (label, a1, a2) for the isogeny classes of abelian surfaces over F_q in
    the LMFDB (through the LMFDB checkout in the configuration).
    """
    db = lmfdb_db(cfg)
    return [(r["label"], int(r["poly"][1]), int(r["poly"][2]))
            for r in db.av_fq_isog.search({"g": 2, "q": q}, ["label", "poly"])]


def kummer_stage(q, datadir, cfg=None):
    """
    kummer_av_classes: abelian surfaces whose wedge square has the same
    transcendental part.  realizable_kummer is -1 when there are none
    (e.g. whenever the geometric Picard rank is below 18) or when the class
    fails the necessary conditions for K3 surfaces (realizable_k3 = -1),
    and 0 otherwise.
    Left null in characteristic 2, when the LMFDB has no abelian surfaces
    over F_q, and when no LMFDB checkout is configured.
    """
    p = ZZ(q).prime_factors()[0]
    if p == 2:
        avs = []
    elif cfg is None or not cfg.lmfdb_path:
        print("q=%s: no LMFDB configured; Kummer columns will be null" % q)
        avs = []
    else:
        avs = lmfdb_abelian_surfaces(q, cfg)
    lookup = defaultdict(list)
    for lab, a1, a2 in avs:
        lookup[tuple(wedge_square_transcendental(q, a1, a2))].append(lab)
    rows = []
    for row in iter_rows(datadir, q, ["label", "geometric_picard_rank", "transcendental_poly",
                                      "realizable_k3"]):
        if p == 2 or not avs:
            rows.append({"label": row["label"], "kummer_av_classes": None, "realizable_kummer": None})
            continue
        labs = sorted(lookup.get(tuple(row["transcendental_poly"]), [])) if row["geometric_picard_rank"] >= 18 else []
        rows.append({"label": row["label"], "kummer_av_classes": labs,
                     "realizable_kummer": 0 if labs and row["realizable_k3"] != -1 else -1})
    os.makedirs(os.path.join(datadir, "kummer"), exist_ok=True)
    write_file(os.path.join(datadir, "kummer", "q%s.txt" % q),
               [("label", "text")] + [(c, ISOG_TYPES[c]) for c in KUMMER_COLUMNS], rows)
    print("q=%s: %s abelian surfaces, %s classes with Kummer candidates"
          % (q, len(avs), sum(1 for r in rows if r["kummer_av_classes"])))


# ----------------------------------------------------------------------
# census: classes known to contain a K3 surface

# Kedlaya and Sutherland, "A census of zeta functions of quartic K3 surfaces
# over F_2": the zeta functions of all smooth quartic surfaces in P^3 over F_2.
QUARTIC_CENSUS_F2 = os.path.join(ROOT, "external", "root-unitary", "k3-quartic-f2", "k3final.txt")


def quartic_census_labels(fname=QUARTIC_CENSUS_F2):
    """
    The labels of the classes over F_2 that contain a smooth quartic surface.
    Each line of the census file is ``[f, [ids]]``, where f(x) = P(x)/(1 - x)
    (a quartic surface has Picard rank at least 1) and the ids number the
    quartics with this zeta function.
    """
    from sage.all import PolynomialRing, ZZ
    Rx = PolynomialRing(ZZ, "x")
    x = Rx.gen()
    labels = set()
    with open(fname) as F:
        for line in F:
            line = line.strip()
            if not line:
                continue
            f = Rx(line[1:line.index(", [")])
            labels.add(poly_to_label((f * (1 - x)).list()))
    return labels


def census_stage(q, datadir, census_file=None):
    """
    Set realizable_k3 = 1 for the classes known to contain a K3 surface over
    F_q.  At present this is the quartic census over F_2, whose classes also
    get polarization degree 4: the hyperplane class of a quartic surface is a
    primitive polarization of degree 4 defined over F_2.  The list is not
    complete (polarization_degrees_complete = False).  The output only has
    rows for those classes; the export keeps the other stages' values (and
    nulls) for the others.  A census class that fails one of the necessary
    conditions of the basic stage would be a bug, so it raises an error.
    """
    if q != 2:
        print("q=%s: no census available" % q)
        return 0
    known = quartic_census_labels(census_file or QUARTIC_CENSUS_F2)
    rows = []
    for row in iter_rows(datadir, q, ["label", "realizable_k3"]):
        if row["label"] in known:
            if row["realizable_k3"] == -1:
                raise RuntimeError("%s contains a smooth quartic but fails a necessary condition" % row["label"])
            rows.append({"label": row["label"], "realizable_k3": 1,
                         "polarization_degrees": [4], "polarization_degrees_complete": False})
    os.makedirs(os.path.join(datadir, "census"), exist_ok=True)
    write_file(os.path.join(datadir, "census", "q%s.txt" % q),
               [("label", "text")] + [(c, ISOG_TYPES[c]) for c in CENSUS_COLUMNS], rows)
    print("q=%s: %s classes in the quartic census, %s of them in the data" % (q, len(known), len(rows)))
    return len(rows)


# ----------------------------------------------------------------------
# refined Artin-Tate conditions (Auel-Shaikh, arXiv:2611.13468)

REFINED_AT_COLUMNS = [("label", "text"), ("realizable_k3", "smallint"), ("refined_at_failure", "integer[]")]


@worker("refinedat")
def refined_at_chunk(cfg, basicfile, outfile):
    from .refined_at import refined_artin_tate_failures
    rows = []
    for row in read_file(basicfile, ["label", "cyclotomic_factors", "transcendental_poly", "realizable_k3"]):
        if row["realizable_k3"] != 0:
            continue
        f = refined_artin_tate_failures(row["cyclotomic_factors"], row["transcendental_poly"])
        if f:
            rows.append({"label": row["label"], "realizable_k3": -1, "refined_at_failure": list(f[0])})
    write_file(outfile, REFINED_AT_COLUMNS, rows)


def refined_at_stage(q, datadir, runner=None):
    """
    The classes that pass the necessary conditions of the basic stage but fail
    a refined Artin-Tate condition (see refined_at.py), with the first failing
    (base change degree l, root of unity order m).  The export sets their
    realizable_k3 (and realizable_kummer, where defined) to -1.  It is an
    error for a class of the census to be among them.
    """
    runner = runner or Runner()
    workdir = os.path.join(datadir, "refined_at", "q%s" % q)
    os.makedirs(workdir, exist_ok=True)
    tasks = []
    outs = []
    for fname in basic_files(datadir, q):
        outfile = os.path.join(workdir, "rat." + os.path.basename(fname))
        outs.append(outfile)
        if not _newer(outfile, [fname]):
            tasks.append((fname, outfile))
    runner.run("refinedat", tasks)
    seen = set()
    rows = []
    for outfile in outs:
        for row in read_file(outfile):
            if row["label"] not in seen:
                seen.add(row["label"])
                rows.append(row)
    census = _load_by_label(os.path.join(datadir, "census", "q%s.txt" % q)) or {}
    clash = [r["label"] for r in rows if r["label"] in census]
    if clash:
        raise RuntimeError("%s classes of the census fail a refined Artin-Tate condition, e.g. %s" % (len(clash), clash[0]))
    rows.sort(key=lambda r: label_sort_key(r["label"]))
    write_file(os.path.join(datadir, "refined_at", "q%s.txt" % q), REFINED_AT_COLUMNS, rows)
    print("q=%s: %s classes ruled out by the refined Artin-Tate conditions" % (q, len(rows)), flush=True)
    return len(rows)


# ----------------------------------------------------------------------
# export

def _load_by_label(fname):
    if not os.path.exists(fname):
        return None
    return {row.pop("label"): row for row in read_file(fname)}


def export_stage(q, datadir):
    fields = None
    fname = os.path.join(datadir, "fields", "q%s.txt" % q)
    if os.path.exists(fname):
        fields = {tuple(r["factor"]): r for r in read_file(fname)}
    else:
        print("q=%s: no fields data; the number field columns will be null" % q)
    extra = {}
    for stage in ("twists", "basechange", "kummer"):
        data = _load_by_label(os.path.join(datadir, stage, "q%s.txt" % q))
        if data is None:
            print("q=%s: no %s data; those columns will be null" % (q, stage))
        else:
            extra[stage] = data
    # only lists the classes ruled out by the refined Artin-Tate conditions
    refined = _load_by_label(os.path.join(datadir, "refined_at", "q%s.txt" % q))
    if refined is None:
        refined = {}
        print("q=%s: no refined Artin-Tate data; realizable_k3 only reflects the basic conditions" % q)
    # only lists the classes known to contain a K3 surface
    census = _load_by_label(os.path.join(datadir, "census", "q%s.txt" % q))
    if census is None:
        census = {}
        if q == 2:
            print("q=2: no census data; realizable_k3 only reflects the necessary conditions")
    outdir = os.path.join(datadir, "export")
    os.makedirs(outdir, exist_ok=True)

    def rows():
        for row in iter_rows(datadir, q):
            facs = [tuple(f) for f, e in row["transcendental_factors"]]
            if fields is not None:
                row["galois_groups"] = [g if g is not None else fields[f]["galois_group"]
                                        for f, g in zip(facs, row["galois_groups"])]
                row["number_field_coeffs"] = [fields[f]["nf_coeffs"] for f in facs]
                row["number_fields"] = [fields[f]["nf_label"] or "" for f in facs]
            for data in extra.values():
                row.update(data[row["label"]])
            if row["label"] in refined:
                if row["label"] in census:
                    raise RuntimeError("%s is in the census but fails a refined Artin-Tate condition" % row["label"])
                row["realizable_k3"] = -1
                if row.get("realizable_kummer") is not None:
                    row["realizable_kummer"] = -1
            row.update(census.get(row["label"], {}))
            yield row
    out = os.path.join(outdir, "k3_fq_isog.q%s.txt" % q)
    write_file(out, ISOG_COLUMNS, rows())
    print("wrote %s" % out)
    return out


@worker("export")
def export_worker(cfg, q, datadir):
    export_stage(int(q), datadir)
