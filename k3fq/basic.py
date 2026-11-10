"""
The basic stage: compute the invariants of K3Poly.compute for every
polynomial in a polynomial file, in chunks run in parallel.

Output: ``basic/q{q}/{stem}.{chunk:05d}.txt`` in the LMFDB upload format,
with columns BASIC_COLUMNS.  Existing chunk files are skipped, so the stage
can be interrupted and resumed.  Chunk files that do not fit the current
polynomial file and chunk size (left by a different enumeration or chunk
size) are an error, since every later stage reads all of basic/q{q}.
"""
import glob
import os
import sys
import time
import traceback

from .columns import BASIC_COLUMNS, ISOG_TYPES
from .invariants import K3Poly
from .parallel import Runner, worker
from .pgformat import count_rows, write_file

CHUNKSIZE = 20000


def chunk_offsets(polyfile, chunksize):
    """
    Yield (index, byte offset, number of lines) for the chunks of ``polyfile``.
    """
    with open(polyfile, "rb") as F:
        index, offset, count, pos = 0, 0, 0, 0
        for line in F:
            if count == chunksize:
                yield index, offset, count
                index, offset, count = index + 1, pos, 0
            count += 1
            pos += len(line)
        if count:
            yield index, offset, count


@worker("basic")
def basic_chunk(cfg, polyfile, offset, count, outfile):
    """
    Compute the basic columns for ``count`` polynomials of ``polyfile``,
    starting at byte ``offset``, and write them to ``outfile``.
    """
    t0 = time.time()
    rows = []
    with open(polyfile, "rb") as F:
        F.seek(int(offset))
        for _ in range(int(count)):
            coeffs = [int(c) for c in F.readline().split()]
            try:
                rows.append(K3Poly(coeffs).compute())
            except Exception:
                sys.stderr.write("Error on %s\n%s" % (coeffs, traceback.format_exc()))
                raise
    write_file(outfile, [(c, ISOG_TYPES[c]) for c in BASIC_COLUMNS], rows)
    print("%s: %s rows in %.1fs" % (os.path.basename(outfile), len(rows), time.time() - t0), flush=True)


def basic_stage(polyfile, datadir, runner=None, chunksize=CHUNKSIZE):
    """
    Run the basic stage on ``polyfile`` (named q{q}_*.txt).  Returns the
    number of polynomials in the file.
    """
    stem = os.path.splitext(os.path.basename(polyfile))[0]
    q = stem.split("_")[0][1:]
    outdir = os.path.join(datadir, "basic", "q%s" % q)
    os.makedirs(outdir, exist_ok=True)
    tasks = []
    total = 0
    expected = set()
    for i, offset, count in chunk_offsets(polyfile, chunksize):
        total += count
        outfile = os.path.join(outdir, "%s.%05d.txt" % (stem, i))
        expected.add(outfile)
        if not os.path.exists(outfile):
            tasks.append((polyfile, offset, count, outfile))
        elif count_rows(outfile) != count:
            raise RuntimeError("%s has %s rows rather than %s: it was made with another chunk size "
                               "or polynomial file; move %s away" % (outfile, count_rows(outfile), count, outdir))
    stale = sorted(set(glob.glob(os.path.join(outdir, "*.txt"))) - expected)
    if stale:
        raise RuntimeError("%s has %s files that do not come from %s (e.g. %s); move them away"
                           % (outdir, len(stale), polyfile, os.path.basename(stale[0])))
    print("%s: %s polynomials, %s chunks to compute" % (os.path.basename(polyfile), total, len(tasks)), flush=True)
    (runner or Runner()).run("basic", tasks)
    return total
