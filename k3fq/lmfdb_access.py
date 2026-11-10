"""
Access to an LMFDB database through a local LMFDB checkout (the [lmfdb]
path in the configuration).  The LMFDB is imported from within the checkout,
as when it is run there: its config.ini and secrets.ini are used, and its log
files (flasklog, slow_queries.log) are written there.
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
    # secrets.ini and the log files are relative to the working directory
    cwd = os.getcwd()
    os.chdir(path)
    try:
        from lmfdb import db
    finally:
        os.chdir(cwd)
    return db
