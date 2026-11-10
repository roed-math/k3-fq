"""
Access to an LMFDB database through a local LMFDB checkout (the [lmfdb]
path in the configuration), which uses that checkout's config.ini.  The
LMFDB's log files (flasklog, slow_queries.log) go to the data directory.
"""
import os
import sys


def lmfdb_db(cfg):
    if cfg is None or not cfg.lmfdb_path:
        raise RuntimeError("this needs an LMFDB checkout: set [lmfdb] path in the configuration")
    path = os.path.abspath(os.path.expanduser(cfg.lmfdb_path))
    if not os.path.isdir(os.path.join(path, "lmfdb")):
        raise RuntimeError("%s is not an LMFDB checkout" % path)
    if path not in sys.path:
        sys.path.insert(0, path)
    cwd = os.getcwd()
    os.makedirs(cfg.datadir, exist_ok=True)
    os.chdir(cfg.datadir)
    try:
        from lmfdb import db
    finally:
        os.chdir(cwd)
    return db
