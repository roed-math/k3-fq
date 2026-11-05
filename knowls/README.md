# Knowls for motives of K3 type over finite fields

These are the knowls used by the `k3fq` section of the LMFDB (`lmfdb/k3fq/` on
the `k3fq` branch, url `/Motive/K3/Fq/`), in the format read by
`db.knowldb.yaml_import`:

- `k3_knowls.yaml`: the definitions.  General K3 notions have ids `k3.*`, and
  notions specific to finite fields have ids `k3.fq.*`.
- `k3_rcs_knowls.yaml`: the source, acknowledgment, citation, completeness and
  reliability pages (`rcs.*.k3.fq`, following the other sections).  The
  acknowledgment and citation pages still have TODOs.

From an LMFDB checkout with write access to the knowls:

```python
from lmfdb.knowledge.knowl import knowldb
knowldb.yaml_import("k3_knowls.yaml", "your_user_id", dryrun=True)   # check first
knowldb.yaml_import("k3_knowls.yaml", "your_user_id")
```

Knowls for interesting classes follow the usual pattern `k3.fq.<label>.top`,
and are listed on `/Motive/K3/Fq/interesting`.

Suggested edit to an existing knowl: the page `/Motive` renders the knowl
`motives.about`, which currently lists only hypergeometric motives.  Add

    <li><a href="{{ url_for('k3fq.index') }}">Motives of K3 type over finite fields</a>

after the hypergeometric line.
