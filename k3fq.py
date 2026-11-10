#!/usr/bin/env sage-python
"""
Build the k3_fq_isog table (motives of K3 type over finite fields) for the
LMFDB.

    sage -python k3fq.py                      run every unfinished stage
    sage -python k3fq.py status               which stages are finished
    sage -python k3fq.py rerun STAGE [Q ...]  redo a stage (and those using it)
    sage -python k3fq.py adopt-markers        keep markers from before 1649994
    sage -python k3fq.py counts Q ...         count classes without listing them

The data directory, the fields and the programs are set in the configuration
file (default: k3fq.toml next to this script); see README.md.
"""
import argparse
import sys


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default=None, help="configuration file (default: k3fq.toml)")
    sub = parser.add_subparsers(dest="command")
    sub.add_parser("run", help="run every unfinished stage (the default)")
    sub.add_parser("status", help="list the stages and whether they are finished")
    p = sub.add_parser("rerun", help="mark a stage and the stages using it as unfinished, then run")
    p.add_argument("stage")
    p.add_argument("q", type=int, nargs="*", help="only for these fields (default: all)")
    sub.add_parser("adopt-markers", help="keep the work recorded by markers from before 1649994 "
                   "(only if the configuration has not changed since)")
    p = sub.add_parser("counts", help="count classes by geometric Picard rank, without listing them")
    p.add_argument("q", nargs="+", help="prime powers, or ranges such as 3-64")
    p.add_argument("-j", "--jobs", type=int, default=None,
                   help="enumeration threads (default: jobs in the configuration)")
    p.add_argument("--max-degree", type=int, default=8,
                   help="largest transcendental degree 2k to count; 22 counts everything (default 8)")
    p.add_argument("--output", default=None, help="also write the markdown table here")
    # used by GNU parallel
    p = sub.add_parser("worker")
    p.add_argument("name")
    p.add_argument("args", nargs="*")

    args = parser.parse_args(argv)
    from k3fq.config import load_config
    cfg = load_config(args.config)
    if args.command == "worker":
        from k3fq.parallel import WORKERS, load_workers
        load_workers()
        WORKERS[args.name](cfg, *args.args)
        return
    if args.command == "counts":
        from k3fq.counts import counts_table, parse_qs
        counts_table(parse_qs(args.q), args.max_degree, threads=args.jobs or cfg.jobs, out=args.output)
        return
    from k3fq.pipeline import Pipeline
    pipeline = Pipeline(cfg)
    if args.command == "status":
        pipeline.status()
        return
    if args.command == "adopt-markers":
        pipeline.adopt_markers()
        return
    if args.command == "rerun":
        pipeline.rerun(args.stage, args.q)
    pipeline.run()


if __name__ == "__main__":
    sys.exit(main())
