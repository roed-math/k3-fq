"""
Reading and writing files in the format used by psycodict's ``copy_from``
and ``reload`` (the LMFDB upload format):

    col1|col2|...
    type1|type2|...
    <blank line>
    val1|val2|...

Arrays use the Postgres text representation ``{1,2,3}``, jsonb is JSON,
booleans are ``t``/``f`` and null is ``\\N``.
"""
import json
import os

SEP = "|"
NULL = r"\N"


def _array_elt(x, base):
    if x is None:
        return "NULL"
    if base == "text":
        s = str(x).replace("\\", "\\\\").replace('"', '\\"')
        return '"%s"' % s
    if base == "boolean":
        return "t" if x else "f"
    if isinstance(x, float):
        return repr(x)
    return str(x)


def _array(x, base):
    if isinstance(x, (list, tuple)):
        return "{" + ",".join(_array(y, base) for y in x) + "}"
    return _array_elt(x, base)


def encode(x, typ):
    if x is None:
        return NULL
    if typ.endswith("[]"):
        return _array(list(x), typ[:-2])
    if typ == "jsonb":
        return json.dumps(x, separators=(",", ":"))
    if typ == "boolean":
        return "t" if x else "f"
    if typ == "text":
        s = str(x)
        if SEP in s or "\n" in s or "\\" in s:
            raise ValueError("unsupported character in %r" % s)
        return s
    if isinstance(x, float):
        return repr(x)
    return str(int(x)) if typ in ("smallint", "integer", "bigint") else str(x)


def _split_array(s):
    """
    Split the inside of a (possibly nested) Postgres array literal.
    """
    parts, depth, cur, quoted, esc = [], 0, [], False, False
    for ch in s:
        if esc:
            cur.append(ch)
            esc = False
        elif ch == "\\":
            esc = True
        elif ch == '"':
            quoted = not quoted
            cur.append(ch)
        elif quoted:
            cur.append(ch)
        elif ch == "{":
            depth += 1
            cur.append(ch)
        elif ch == "}":
            depth -= 1
            cur.append(ch)
        elif ch == "," and depth == 0:
            parts.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    if cur or parts:
        parts.append("".join(cur))
    return parts


def _decode_elt(s, base):
    if s == "NULL":
        return None
    if s.startswith("{"):
        return [_decode_elt(t, base) for t in _split_array(s[1:-1])]
    if base == "text":
        if s.startswith('"'):
            s = s[1:-1]
        return s.replace('\\"', '"').replace("\\\\", "\\")
    if base == "boolean":
        return s == "t"
    if base in ("float8", "real", "double precision"):
        return float(s)
    return int(s)


def decode(s, typ):
    if s == NULL:
        return None
    if typ.endswith("[]"):
        return _decode_elt(s, typ[:-2])
    if typ == "jsonb":
        return json.loads(s)
    if typ == "boolean":
        return s == "t"
    if typ == "text":
        return s
    if typ in ("float8", "real", "double precision"):
        return float(s)
    return int(s)


def write_file(filename, columns, rows):
    """
    Write ``rows`` (dictionaries) to ``filename``, atomically.  ``columns``
    is a list of (name, type) pairs.
    """
    tmp = filename + ".tmp"
    with open(tmp, "w") as F:
        F.write(SEP.join(c for c, _ in columns) + "\n")
        F.write(SEP.join(t for _, t in columns) + "\n\n")
        for row in rows:
            F.write(SEP.join(encode(row.get(c), t) for c, t in columns) + "\n")
    os.replace(tmp, filename)


def read_file(filename):
    """
    Iterate over the rows of ``filename`` as dictionaries.
    """
    with open(filename) as F:
        cols = F.readline().rstrip("\n").split(SEP)
        types = F.readline().rstrip("\n").split(SEP)
        F.readline()
        for line in F:
            vals = line.rstrip("\n").split(SEP)
            yield {c: decode(v, t) for c, t, v in zip(cols, types, vals)}


def read_header(filename):
    with open(filename) as F:
        cols = F.readline().rstrip("\n").split(SEP)
        types = F.readline().rstrip("\n").split(SEP)
    return list(zip(cols, types))
