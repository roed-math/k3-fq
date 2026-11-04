"""
The column list in k3fq/columns.py must match schemas/k3_fq_isog.md.
"""
import os
import re

from k3fq.columns import ISOG_COLUMNS

SCHEMA = os.path.join(os.path.dirname(os.path.dirname(__file__)), "schemas", "k3_fq_isog.md")


def schema_columns():
    cols = []
    with open(SCHEMA) as F:
        for line in F:
            m = re.match(r"^\| ([a-z0-9_]+) \| ([a-z0-9 \[\]]+) \|", line)
            if m and m.group(1) != "Column":
                cols.append((m.group(1), m.group(2)))
    return cols


def test_columns_match_schema():
    assert schema_columns() == ISOG_COLUMNS
