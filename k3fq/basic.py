"""
The basic stage: compute the invariants of K3Poly.compute for every
polynomial in a polynomial file, in parallel over chunks.

Output: ``basic/q{q}/{stem}.{chunk:05d}.txt`` in the LMFDB upload format,
with columns BASIC_COLUMNS.  Existing chunk files are skipped, so the stage
can be interrupted and resumed.
"""
import os
import sys
import time
import traceback
from multiprocessing import Pool

from .columns import BASIC_COLUMNS, ISOG_TYPES
from .invariants import K3Poly
from .pgformat import write_file

CHUNKSIZE = 5000


def _chunks(polyfile, chunksize):
    chunk = []
    with open(polyfile) as F:
        for line in F:
            chunk.append(line)
            if len(chunk) == chunksize:
                yield chunk
                chunk = []
    if chunk:
        yield chunk


def _process(args):
    lines, outfile = args
    t0 = time.time()
    rows = []
    for line in lines:
        coeffs = [int(c) for c in line.split()]
        try:
            rows.append(K3Poly(coeffs).compute())
        except Exception:
            sys.stderr.write("Error on %s\n%s" % (coeffs, traceback.format_exc()))
            raise
    write_file(outfile, [(c, ISOG_TYPES[c]) for c in BASIC_COLUMNS], rows)
    return outfile, len(rows), time.time() - t0


def basic_stage(polyfile, datadir, jobs=8, chunksize=CHUNKSIZE):
    stem = os.path.splitext(os.path.basename(polyfile))[0]
    q = stem.split("_")[0][1:]
    outdir = os.path.join(datadir, "basic", "q%s" % q)
    os.makedirs(outdir, exist_ok=True)

    def tasks():
        for i, chunk in enumerate(_chunks(polyfile, chunksize)):
            outfile = os.path.join(outdir, "%s.%05d.txt" % (stem, i))
            if not os.path.exists(outfile):
                yield chunk, outfile

    t0 = time.time()
    total = 0
    with Pool(jobs) as pool:
        for outfile, n, t in pool.imap_unordered(_process, tasks()):
            total += n
            print("%s: %s rows in %.1fs (total %s, %.1fs elapsed)"
                  % (os.path.basename(outfile), n, t, total, time.time() - t0), flush=True)
    return total
