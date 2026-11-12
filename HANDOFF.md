# Handoff: motives of K3 type over finite fields (LMFDB)

Written 2026-11-05 by the claudeF agent on mimir, for whoever picks this up next.  The project owner is David Roe.  Shared copy: `/home/agents/fleet-shared/k3-fq/HANDOFF.md`; repository copy: `HANDOFF.md` in roed-math/k3-fq.

## Goal

This is a new LMFDB section, "Motives of K3 type" over F_q.  Each row is a degree-22 polynomial P(T) in Z[T] with constant term q, leading coefficient ±q and all roots on the unit circle.  It is modeled on abelian varieties over F_q, and grouped under Motives in the sidebar.  David's original notes are in `/home/claudeF/k3_fq.md` on mimir.

## Where everything is

| Piece | Location | State |
| --- | --- | --- |
| Compute repo | https://github.com/roed-math/k3-fq (public), branch `main`, head `2ed87aa`; mimir clone `/scratch/claudeF/k3-fq` | pushed, clean |
| Schemas | `schemas/k3_fq_isog.md`, `schemas/k3_fq_lattices.md` in the compute repo | identical to PR #75 |
| Schema PR | https://github.com/jvoight/k3s-lmfdb/pull/75 (from roed-math branch `k3-fq-schemas`, head `d004a16`; mimir worktree `/scratch/claudeF/k3s-lmfdb-schemas`) | open, no review yet |
| LMFDB front end | roed-math/lmfdb branch `k3fq`, head `78bb194a3`, 9 commits on upstream main `8aedadac9`; mimir worktree `/scratch/claudeF/lmfdb-k3fq` | pushed; **no PR to LMFDB opened** |
| Knowls | `knowls/k3_knowls.yaml` (27 knowls: everything the section links to, `k3.*` for general K3 notions and `k3.fq.*` for finite fields) and `knowls/k3_rcs_knowls.yaml` (the five `rcs.*.k3.fq` pages; acknowledgment and citation still have TODOs), both for `db.knowldb.yaml_import`; see `knowls/README.md` | not in any database yet |
| Test upload | `test_upload/k3_fq_isog_test.txt` (13 rows) and `test_upload/create_k3_fq_isog.py` (the `create_table` call, generated from the schema by `scripts/make_create_table.py`) | never run against a real database |
| Kedlaya's root-unitary | git submodule `external/root-unitary` (pinned at `d046f6c`); `make` builds `bin/weilpoly` from `c-script/weilpoly.c` against Sage's FLINT | |
| Sample data (mimir) | `/scratch/claudeF/work/k3fq_sample_data/` (pipeline dirs), combined export `/scratch/claudeF/work/k3fq_sample/k3_fq_isog.txt` (4,809 rows, q = 2, 3, 4) | readable by all; only claudeF can write.  The same sample with the new labels is `/scratch/claudeM/work/k3fq_sample/k3_fq_isog.txt` (pipeline dirs in `/scratch/claudeM/work/k3fq_sample_data/`) |
| Front-end harness (mimir) | `/scratch/claudeF/work/k3fq_frontend/`: `run_harness.py [--sample=FILE]` (Flask test client against an in-memory fake table), `run_unit_tests.py --sample=FILE`, `serve.py` (port 37791) | run from the LMFDB worktree with `sage -python` |

## Pipeline (compute repo)

`sage -python k3fq.py` runs the whole pipeline (`enumerate`, `basic`, `fields`, `twists`, `kummer`, `census`, `basechange`, `export`, `combine`) as set up by the configuration file `k3fq.toml`: data directory, jobs, programs, LMFDB checkout, and per field q the rank cutoff and expected class count.  Finished stages leave markers in `DATA/done/`, so rerunning the command resumes; `status` and `rerun STAGE [Q...]` manage them.  Parallel work goes through GNU parallel, one `sage -python k3fq.py --config FILE worker NAME ARGS` process per task (`k3fq/parallel.py`).  `counts` is separate.  The README has the details.  Tests: `make test` (93 pass; a fresh `git clone --recursive` passes too).

Key design points:
- Labels are `q.ρ.c.t` (`k3fq/labels.py`): the geometric Picard rank, the index `c` of the cyclotomic part among those of degree ρ̄, and the coefficients of the transcendental part in balanced base 26 (`n` = 0), all padded to one width and concatenated.  They replaced `q.s.c1_..._c11` on 2026-11-05, because large cyclotomic parts made those 25–58 characters long (median 32 in the sample; the new median is 11).  The geometric class representative of a complete twist class always has `c = a`.
- Sign −1 polynomials are (1 − T²) times the sign +1 polynomials of degree 20.
- The enumeration multiplies transcendental parts of degree ≤ 22 − R by every cyclotomic product of the complementary degree (R = `min_geometric_rank`; R = 0 for all of F_2).
- Galois groups are certified as C2 ≀ Gal(h) via Frobenius cycle types, with PARI's `polgalois` for the real polynomial h; Magma is the fallback.
- Twist classes are found by hashing the multiset {Mθ} with M = lcm{m : φ(m) ≤ 44}.  Base change uses power sums.
- Kummer candidates come from `av_fq_isog`, and number-field labels from `nf_fields` (looked up by the `polredabs` polynomial), both through the LMFDB checkout in the configuration (`~roed/lmfdb` on lovelace, which connects locally).
- `counts` counts classes per geometric Picard rank from the driver's totals, without listing them.  This uses N(q, ρ̄ ≥ R) = Σ T(q,2k)·C(22−2k), with C(22) = 60,280.

Not implemented (exported as null):
- `polarization_degrees` beyond the census, which gives degree 4 for its classes, with `polarization_degrees_complete` false.  Taelman's theorem would give complete lists for ordinary classes; results that need potential semistable reduction (p ≤ 3 in our data) are to be kept out of this column (decided 2026-11-06).
- The `k3_fq_lattices` computation.  A prototype in `lattices/` (Julia + Oscar, following Auel–Shaikh arXiv:2611.13468: rootless coinvariant lattices with discriminants bounded by their Theorem 4.1, glued to the invariant lattice, then Artin–Tate, Brauer-map and ℓ-adic checks) was tested on a sample of q = 2, ρ̄ = 4 (2026-11-12); see `lattices/README.md`.
- The `refined_at` stage (2026-11-12) adds the refined Artin–Tate conditions of Auel–Shaikh (arXiv:2611.13468) with the bounds of Auel–Engel–Kedlaya–Shaikh (arXiv:2611.08356) as a necessary condition for `realizable_k3`.  No class of the census fails it; on a sample over F_2 with ρ̄ = 4 it rules out about a third of the classes that pass the basic conditions.  `k3fq.py update-file` writes an update file for a table loaded before it existed.
- Positive realizability beyond the quartic census over F_2 (the `census` stage sets `realizable_k3 = 1` for its 52,755 classes; all of them pass our necessary conditions).  Otherwise only necessary conditions are used, giving −1.

## Open decisions and next steps

1. **Cutoffs (decided 2026-11-05; q = 9, 11 added 2026-11-12, for the pointless K3 surfaces).**  All classes for q = 2; geometric Picard rank ρ̄ ≥ 16 for q = 3, ≥ 18 for q = 4, and ≥ 20 for q = 5, 7, 8, 9, 11 (recorded in `rcs.cande.k3.fq`; 15,665,490 classes, with 1,037,304 for q = 9 and 1,281,560 for q = 11).  The q ≥ 3 part is 7,868,568 rows (about 1.44 times F_2), so `k3_fq_isog` has about 13.35M rows in all; all of F_2 is about 5.48M classes (5,403,212 with ρ̄ ≥ 2, plus a ρ̄ = 0 part of about 80k, not yet pinned down).  The sign −1 count over F_2 is 1,485,591.  Cutoffs must be nondecreasing in q; the primitivity logic relies on this.  Exact counts for all q ≤ 64 are reproduced by `k3fq.py counts 3-64 --max-degree 6`.
2. **Lovelace run (k3_fq_isog only).**  David runs it; claudeF/claudeM cannot SSH there.  `git clone --recursive`, `make`, `make test`; check `k3fq.toml` (data `/scratch/k3`, `jobs = 0` for one per core, Magma, `polredabs`, LMFDB checkout `~roed/lmfdb`); then `nohup sage -python k3fq.py > k3fq.log 2>&1 &` and `sage -python k3fq.py status` to follow it.  The output is `/scratch/k3/k3_fq_isog.txt`.  The lattice data is independent and can be added later with `update_from_file`.  mimir is shared and often heavily loaded (load 100–160 seen), so keep jobs there tiny and niced.
3. **Database test.**  Run `create_k3_fq_isog.py`, then `copy_from` the 13-row file, then `K3FqStats().setup()`.  Check real Postgres behaviour of `$startswith` on `poly_str`, `$overlaps` on `galois_groups`, and null ordering.  Add indexes on (q, poly) and on `poly_str` with `text_pattern_ops`.
4. **LMFDB integration.**
   - Import the knowls with `db.knowldb.yaml_import` (see `knowls/README.md`).  `rcs.ack.k3.fq` and `rcs.cite.k3.fq` still have TODOs.
   - Add the k3fq link to the `motives.about` knowl (it renders `/Motive`).
   - Open the LMFDB PR when David says so.
5. **Mathematical conventions to confirm with David.**
   - "Igusa–Makai" in the notes was read as Gushel–Mukai.
   - The Verra fourfold point count assumes the 22-dimensional piece contains one of the three classes pulled back from P²×P².
   - Both are flagged in the schema notes.
6. **Performance at scale (measured 2026-11-09 on mimir at load 60–90, so pessimistic).**  Cutoffs as in item 1: 13,351,780 classes.  CPU-hours by stage:

   | stage | cost | total |
   | --- | --- | --- |
   | enumerate | F_2 degrees 18, 20, 22: 0.24, 2.0, 5.6 CPU-h (528,687, 1,485,591 and 3,992,467 polynomials, so F_2 has exactly 5,478,058 classes, 74,846 of them with ρ̄ = 0); q ≥ 3 minutes | ~8 |
   | basic | 1.8 ms per class when its transcendental part has been seen before (cached), plus ~30 ms per distinct transcendental part (~154k over F_2) | ~5 |
   | fields | `polredbest` and `polredabs` ≤ 0.05 s per factor even in degree 22 (discriminants have ≤ ~40 digits); Magma for the 30–40% of factors of degree 12, 16, 20 that the certificate misses (0 in degrees 18, 22), 0.2 s each, ~17k factors | ~3 |
   | twists: keys | 0.12 ms per class (cached by transcendental part) | ~1 |
   | twists: pairwise | only for ρ̄ ≤ 4 over F_2: every transcendental part pairs with its quadratic twist, so twist classes have 2 × C(ρ̄) members, which is ≤ 50 only for ρ̄ ≤ 4.  The minimal twist degrees are computed separately for the cyclotomic parts (exact, cached by pair) and the transcendental parts (numerical, verified by exact base change of the transcendental parts, cached by pair), and base-change labels are assembled from the two; ~0.01 s per class at ρ̄ = 4 (441k classes), less below | ~2 |
   | basechange, kummer, census | streamed; base change only for subfield classes of large enough rank | < 1 |
   | export | 0.3 ms per row, ~1.5 KB per row | ~2, ~20 GB |

   So ~20 CPU-h in all for `k3_fq_isog`, a few hours on one large VM.  The enumeration writes the classes grouped by transcendental part, so the caches in `basic` and `twists` hit almost every time.  The export holds the F_2 stage data in memory (several GB, so `export_jobs = 1` by default); 64 GB RAM and ~60 GB of disk are enough (the basic stage, the per-field exports and the combined `k3_fq_isog.txt` are about 1.1 KB per class each).  Magma (for the `fields` stage) may not be licensed on cloud VMs.

   **Lattices** (pointless classes imported from Simon Brandhorst, E ≤ 6 computed) are the unknown: 1.9M classes with E ≤ 6 and about 55M candidate genera after two necessary conditions, while Brandhorst's pointless run averaged ~0.9 CPU-h per genus (100 CPU-years for ~10^6 genera).  Run a pilot with his code on a few hundred E ≤ 2 and E ≤ 6 classes before committing.

## Reviews already addressed

David passed along two GPT reviews, and both are fully addressed.

First review:
- Untracked support files: the `.gitignore` pattern was unanchored.
- `is_ordinary` is now p ∤ a₁.
- Twist labels are deduplicated across overlapping inputs.
- Counts ranges and characteristics in the front end.
- Nested delimiters in number-field names.
- The lattice schema is now K3-only; fourfold lattice data needs separate tables.
- Fourfold motives are H⁴(2), weight 0.
- Taelman's existence hypotheses are stated.
- The primitivity item was withdrawn, given nondecreasing cutoffs.

Second review:
- Counts accepts arithmetic expressions by reusing `parse_ints`.
- Knowl wording.
- Hassett's "admissible discriminants".
- A weight-0 code comment.

## Gotchas

- In Sage, `QQ(a/b) % 1` is modular inversion, not the fractional part; use `x - x.floor()`.
- psycodict compares upload-file header types with Postgres `udt_name::regtype` names: write `double precision[]`, not `float8[]`.
- There is no Postgres on mimir.  Running binaries extracted from a downloaded pip wheel was blocked by the permission classifier; don't retry that.  Use the fake-table harness, or ask David for a database.
- Firefox (a snap) can only write screenshots into `$HOME`.
- `sage -python` has flask and psycodict, so use it for the LMFDB front end.  The `/scratch/claudeF/sagevenv` venv has no Sage.
- The `@count_wrap` and `@search_wrap` decorators in `lmfdb/k3fq/main.py` apply to the next `def`.  Don't insert helpers between a decorator and its function.
- Ask David before pushing or opening PRs.  He approved everything pushed so far.
