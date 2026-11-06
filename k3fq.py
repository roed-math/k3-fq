#!/usr/bin/env sage-python
"""
Build tables of motives of K3 type over finite fields for the LMFDB.

Run with ``sage -python k3fq.py <stage> ...``; see README.md.
"""
import argparse
import sys


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", default="data", help="data directory (default: data)")
    parser.add_argument("-j", "--jobs", type=int, default=8, help="parallel processes")
    sub = parser.add_subparsers(dest="stage", required=True)

    p = sub.add_parser("enumerate", help="enumerate K3-type polynomials with root-unitary")
    p.add_argument("q", type=int)
    p.add_argument("--min-geometric-rank", type=int, default=None,
                   help="only polynomials with geometric Picard rank at least this "
                        "(enumerates transcendental parts of small degree)")

    p = sub.add_parser("basic", help="compute the invariants of each polynomial")
    p.add_argument("polyfiles", nargs="+", help="files written by the enumerate stage")
    p.add_argument("--chunksize", type=int, default=5000)

    for name, helptext in [("fields", "Galois groups missed by the certificate (Magma) and reduced field polynomials"),
                           ("twists", "twist classes and minimal twist degrees"),
                           ("kummer", "Kummer candidates from abelian surfaces in the LMFDB"),
                           ("census", "classes known to contain a K3 surface (the quartic census over F_2)"),
                           ("export", "assemble the LMFDB upload file")]:
        p = sub.add_parser(name, help=helptext)
        p.add_argument("q", type=int, nargs="+")
    sub.add_parser("basechange", help="primitive models (uses every q in the data directory)")
    p = sub.add_parser("counts", help="count classes by geometric Picard rank, without listing them")
    p.add_argument("q", nargs="+", help="prime powers, or ranges such as 3-64")
    p.add_argument("--max-degree", type=int, default=8,
                   help="largest transcendental degree 2k to count; 22 counts everything (default 8)")
    p.add_argument("--output", default=None, help="also write the markdown table here")

    args = parser.parse_args(argv)
    if args.stage == "enumerate":
        from k3fq.enumerate import enumerate_full, enumerate_elevated
        import os
        outdir = os.path.join(args.data, "polys")
        if args.min_geometric_rank is None:
            enumerate_full(args.q, outdir, threads=args.jobs)
        else:
            enumerate_elevated(args.q, args.min_geometric_rank, outdir, threads=args.jobs)
    elif args.stage == "basic":
        from k3fq.basic import basic_stage
        for f in args.polyfiles:
            basic_stage(f, args.data, jobs=args.jobs, chunksize=args.chunksize)
    elif args.stage == "counts":
        from k3fq.counts import counts_table, parse_qs
        counts_table(parse_qs(args.q), args.max_degree, threads=args.jobs, out=args.output)
    elif args.stage == "basechange":
        from k3fq.stages import basechange_stage
        basechange_stage(args.data)
    else:
        from k3fq import stages
        func = getattr(stages, args.stage + "_stage")
        for q in args.q:
            if args.stage in ("fields", "twists"):
                func(q, args.data, jobs=args.jobs)
            else:
                func(q, args.data)


if __name__ == "__main__":
    sys.exit(main())
