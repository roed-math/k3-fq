"""
The stages after ``basic``, and the final export.  All intermediate files
live under the data directory:

    polys/q{q}_*.txt          output of ``enumerate`` (one polynomial per line)
    basic/q{q}/*.txt          output of ``basic`` (BASIC_COLUMNS)
    fields/q{q}.txt           Galois groups and number fields, keyed by factor
    twists/q{q}/keys.*.txt    geometric keys (intermediate)
    twists/q{q}.txt           TWIST_COLUMNS
    basechange/q{q}.txt       BASECHANGE_COLUMNS
    kummer/q{q}.txt           KUMMER_COLUMNS
    census/q{q}.txt           CENSUS_COLUMNS, for the classes known to contain a K3 surface
    export/k3_fq_isog.q{q}.txt  all columns, ready for psycodict's copy_from
"""
import glob
import hashlib
import json
import os
import subprocess
import tempfile
import time
from collections import defaultdict
from multiprocessing import Pool

from sage.all import ZZ, pari

from .columns import (ISOG_COLUMNS, ISOG_TYPES, TWIST_COLUMNS, BASECHANGE_COLUMNS,
                      KUMMER_COLUMNS, CENSUS_COLUMNS)
from .cyclotomic import R, split_cyclotomic
from .labels import label_sort_key, poly_to_label
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
# fields: Galois groups that the certificate missed (Magma) and reduced
# defining polynomials (PARI polredbest)

def _polredbest(f):
    g = pari(R(f)).polredbest()
    return [int(c) for c in pari.Vecrev(g)]


def fields_stage(q, datadir, jobs=8, magma="magma"):
    factors = set()
    missing = set()
    for row in iter_rows(datadir, q, ["transcendental_factors", "galois_groups"]):
        for (f, e), g in zip(row["transcendental_factors"], row["galois_groups"]):
            factors.add(tuple(f))
            if g is None:
                missing.add(tuple(f))
    print("q=%s: %s distinct transcendental factors, %s without a Galois group"
          % (q, len(factors), len(missing)))
    galois = {}
    if missing:
        galois = magma_galois_groups(sorted(missing), magma)
    factors = sorted(factors)
    with Pool(jobs) as pool:
        reduced = pool.map(_polredbest, factors, chunksize=200)
    out = os.path.join(datadir, "fields", "q%s.txt" % q)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    rows = [{"factor": list(f), "galois_group": galois.get(f), "nf_label": None,
             "nf_coeffs": g} for f, g in zip(factors, reduced)]
    write_file(out, [("factor", "jsonb"), ("galois_group", "text"), ("nf_label", "text"),
                     ("nf_coeffs", "jsonb")], rows)
    return len(rows)


def magma_galois_groups(polys, magma="magma"):
    """
    Transitive labels of the Galois groups of ``polys`` (coefficient tuples),
    computed by magma/galois_groups.m.
    """
    script = os.path.join(ROOT, "magma", "galois_groups.m")
    with tempfile.TemporaryDirectory() as tmp:
        infile = os.path.join(tmp, "in.txt")
        outfile = os.path.join(tmp, "out.txt")
        with open(infile, "w") as F:
            for f in polys:
                F.write(",".join(map(str, f)) + "\n")
        t0 = time.time()
        subprocess.run([magma, "-b", "infile:=%s" % infile, "outfile:=%s" % outfile, script],
                       check=True, stdout=subprocess.DEVNULL)
        print("Magma: %s Galois groups in %.1fs" % (len(polys), time.time() - t0))
        ans = {}
        with open(outfile) as F:
            for line in F:
                coeffs, label = line.split()
                ans[tuple(int(c) for c in coeffs.split(","))] = label
    assert len(ans) == len(polys)
    return ans


# ----------------------------------------------------------------------
# twists

_TWIST_INPUT = ["label", "q", "cyclotomic_factors", "geometric_picard_rank", "transcendental_poly",
                "transcendental_factors"]


def _keys_for_chunk(args):
    fname, outfile = args
    rows = []
    for row in read_file(fname):
        key = tw.geometric_key(row)
        rows.append({"label": row["label"],
                     "key": hashlib.sha1(key.encode()).hexdigest()})
    write_file(outfile, [("label", "text"), ("key", "text")], rows)
    return len(rows)


def _class_twists(members):
    return tw.twist_columns(members, tw.class_twists(members))


def twists_stage(q, datadir, jobs=8):
    keydir = os.path.join(datadir, "twists", "q%s" % q)
    os.makedirs(keydir, exist_ok=True)
    tasks = []
    for fname in basic_files(datadir, q):
        outfile = os.path.join(keydir, "keys." + os.path.basename(fname))
        if not os.path.exists(outfile):
            tasks.append((fname, outfile))
    t0 = time.time()
    with Pool(jobs) as pool:
        for n in pool.imap_unordered(_keys_for_chunk, tasks):
            pass
    print("q=%s: geometric keys computed (%.1fs)" % (q, time.time() - t0))
    # sets: overlapping enumeration files (e.g. full and elevated rank) share labels
    classes = defaultdict(set)
    for fname in sorted(glob.glob(os.path.join(keydir, "keys.*.txt"))):
        for row in read_file(fname):
            classes[row["key"]].add(row["label"])
    small = {lab for labs in classes.values() if 1 < len(labs) <= tw.MAX_TWISTS for lab in labs}
    data = {r["label"]: r for r in iter_rows(datadir, q, _TWIST_INPUT) if r["label"] in small}
    out = {}
    work = []
    for labs in classes.values():
        labs = sorted(labs, key=label_sort_key)
        if 1 < len(labs) <= tw.MAX_TWISTS:
            work.append([data[lab] for lab in labs])
        else:
            out.update(tw.twist_columns([{"label": lab} for lab in labs],
                                        None if len(labs) > 1 else {labs[0]: []}))
    with Pool(jobs) as pool:
        for res in pool.imap_unordered(_class_twists, work):
            out.update(res)
    rows = [dict(label=lab, **cols) for lab, cols in sorted(out.items(), key=lambda kv: label_sort_key(kv[0]))]
    write_file(os.path.join(datadir, "twists", "q%s.txt" % q),
               [("label", "text")] + [(c, ISOG_TYPES[c]) for c in TWIST_COLUMNS], rows)
    print("q=%s: %s twist classes (%.1fs)" % (q, len(classes), time.time() - t0))
    return len(rows)


# ----------------------------------------------------------------------
# primitive models (needs the basic stage for all relevant q)

def basechange_stage(datadir):
    """
    Primitive models.  Base change preserves the geometric Picard rank, so a
    class over a subfield F_{q0} can only be a model of a class in the data
    over F_q if its rank is at least the smallest rank over F_q; only those
    are base changed.  (Since the cutoffs are nondecreasing in q, the data
    over F_{q0} contains all of them.)  The classes over F_q are streamed.
    """
    qs = available_qs(datadir)
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


def lmfdb_abelian_surfaces(q):
    """
    (label, a1, a2) for the isogeny classes of abelian surfaces over F_q in
    the LMFDB (read from the public devmirror).
    """
    import psycopg
    conn = psycopg.connect(host="devmirror.lmfdb.xyz", user="lmfdb", password="lmfdb",
                           dbname="lmfdb")
    with conn.cursor() as cur:
        cur.execute("SELECT label, poly FROM av_fq_isog WHERE g = 2 AND q = %s", (q,))
        ans = [(lab, int(poly[1]), int(poly[2])) for lab, poly in cur.fetchall()]
    conn.close()
    return ans


def kummer_stage(q, datadir):
    """
    kummer_av_classes: abelian surfaces whose wedge square has the same
    transcendental part.  realizable_kummer is -1 when there are none
    (e.g. whenever the geometric Picard rank is below 18) or when the class
    fails the necessary conditions for K3 surfaces (realizable_k3 = -1),
    and 0 otherwise.
    Left null in characteristic 2 and when the LMFDB has no abelian surfaces
    over F_q.
    """
    p = ZZ(q).prime_factors()[0]
    avs = lmfdb_abelian_surfaces(q) if p != 2 else []
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
# export

def _load_by_label(fname):
    if not os.path.exists(fname):
        return None
    return {row.pop("label"): row for row in read_file(fname)}


def export_stage(q, datadir):
    fields = {}
    fname = os.path.join(datadir, "fields", "q%s.txt" % q)
    if os.path.exists(fname):
        fields = {tuple(r["factor"]): r for r in read_file(fname)}
    extra = {}
    for stage in ("twists", "basechange", "kummer"):
        data = _load_by_label(os.path.join(datadir, stage, "q%s.txt" % q))
        if data is None:
            print("q=%s: no %s data; those columns will be null" % (q, stage))
        else:
            extra[stage] = data
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
            if fields:
                row["galois_groups"] = [g if g is not None else fields[f]["galois_group"]
                                        for f, g in zip(facs, row["galois_groups"])]
                row["number_field_coeffs"] = [fields[f]["nf_coeffs"] for f in facs]
                row["number_fields"] = [fields[f]["nf_label"] or "" for f in facs]
            for data in extra.values():
                row.update(data[row["label"]])
            row.update(census.get(row["label"], {}))
            yield row
    out = os.path.join(outdir, "k3_fq_isog.q%s.txt" % q)
    write_file(out, ISOG_COLUMNS, rows())
    print("wrote %s" % out)
    return out
