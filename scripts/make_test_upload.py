"""
Make a small test file for uploading k3_fq_isog to an LMFDB database.

It takes a handful of varied rows from an existing export and adds two
supersingular classes over F_2 that the LMFDB tests (lmfdb/k3fq/test_k3fq.py)
look for.  The later-stage columns for those two are filled in by hand: they
are primitive (F_2 is a prime field), and their twist classes are left as
not computed.

    sage -python scripts/make_test_upload.py EXPORT_FILE OUTPUT_FILE
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from k3fq.columns import ISOG_COLUMNS  # noqa: E402
from k3fq.cyclotomic import R, T  # noqa: E402
from k3fq.invariants import K3Poly  # noqa: E402
from k3fq.pgformat import read_file, write_file  # noqa: E402


def supersingular_rows():
    rows = []
    for P in [2 * (1 - T) ** 22, 2 * (1 - T) ** 21 * (1 + T)]:
        row = K3Poly(R(P).list()).compute()
        row.update(number_fields=[], number_field_coeffs=[], galois_groups=[],
                   is_primitive=True, primitive_models=None,
                   geometric_class=None, twist_count=None, twists=None, max_twist_degree=None,
                   kummer_av_classes=None, realizable_kummer=None)
        rows.append(row)
    return rows


def pick(rows):
    """
    A dozen or so rows with different features, in a fixed order.
    """
    by_label = {r["label"]: r for r in rows}
    chosen = []

    def take(pred):
        for r in rows:
            if r["label"] not in {c["label"] for c in chosen} and pred(r):
                chosen.append(r)
                return r

    generic = take(lambda r: r["q"] == 2 and r["geometric_picard_rank"] <= 6 and r["twists"])
    if generic:
        # one of its twists, if present
        for lab, bc, n in generic["twists"]:
            if lab in by_label:
                chosen.append(by_label[lab])
                break
    take(lambda r: r["q"] == 2 and r["sign"] == -1 and r["realizable_k3"] == 0)
    take(lambda r: r["q"] == 2 and r["sign"] == 1 and r["is_ordinary"] and r["realizable_k3"] == 0)
    take(lambda r: r["q"] == 2 and r["realizable_k3"] == 1)   # in the quartic census
    take(lambda r: r["q"] == 2 and r["is_pointless"])
    take(lambda r: r["q"] == 2 and (r["height"] or 0) >= 3)
    take(lambda r: r["q"] == 3 and r["kummer_av_classes"])
    take(lambda r: r["q"] == 3 and r["sign"] == -1)
    nonprim = take(lambda r: r["q"] == 4 and r["primitive_models"])
    if nonprim:
        for lab in nonprim["primitive_models"]:
            if lab in by_label and by_label[lab] not in chosen:
                chosen.append(by_label[lab])
    take(lambda r: r["q"] == 4 and not r["has_k3_newton_polygon"])
    take(lambda r: r["q"] == 4 and r["transcendental_irreducible"] is False)
    return chosen


def main(infile, outfile):
    rows = list(read_file(infile))
    out = supersingular_rows() + pick(rows)
    write_file(outfile, ISOG_COLUMNS, out)
    for r in out:
        print(r["label"])
    print("%s rows written to %s" % (len(out), outfile))


if __name__ == "__main__":
    main(*sys.argv[1:])
