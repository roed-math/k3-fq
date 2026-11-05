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
| Knowl drafts (37) | `knowls/k3_fq_knowls.md` in the compute repo | not in any database yet |
| Test upload | `test_upload/k3_fq_isog_test.txt` (12 rows) and `test_upload/create_k3_fq_isog.py` (the `create_table` call, generated from the schema by `scripts/make_create_table.py`) | never run against a real database |
| Kedlaya's root-unitary | git submodule `external/root-unitary` (pinned at `d046f6c`); `make` builds `bin/weilpoly` from `c-script/weilpoly.c` against Sage's FLINT | |
| Sample data (mimir) | `/scratch/claudeF/work/k3fq_sample_data/` (pipeline dirs), combined export `/scratch/claudeF/work/k3fq_sample/k3_fq_isog.txt` (4,809 rows, q = 2, 3, 4) | readable by all; only claudeF can write |
| Front-end harness (mimir) | `/scratch/claudeF/work/k3fq_frontend/`: `run_harness.py [--sample=FILE]` (Flask test client against an in-memory fake table), `run_unit_tests.py --sample=FILE`, `serve.py` (port 37791) | run from the LMFDB worktree with `sage -python` |

## Pipeline (compute repo)

The stages are `enumerate`, `basic`, `fields`, `twists`, `basechange`, `kummer`, `export` and `counts`; run them as `sage -python k3fq.py ...`.  The README has the details.  Tests: `make test` (55 pass; a fresh `git clone --recursive` passes too).

Key design points:
- Sign −1 polynomials are (1 − T²) times the sign +1 polynomials of degree 20.
- `--min-geometric-rank R` multiplies transcendental parts of degree ≤ 22 − R by every cyclotomic product of the complementary degree.
- Galois groups are certified as C2 ≀ Gal(h) via Frobenius cycle types, with PARI's `polgalois` for the real polynomial h; Magma is the fallback.
- Twist classes are found by hashing the multiset {Mθ} with M = lcm{m : φ(m) ≤ 44}.  Base change uses power sums.
- Kummer candidates come from `av_fq_isog` on the public devmirror.
- `counts` counts classes per geometric Picard rank from the driver's totals, without listing them.  This uses N(q, ρ̄ ≥ R) = Σ T(q,2k)·C(22−2k), with C(22) = 60,280.

Not implemented (exported as null):
- LMFDB number-field labels (needs polredabs and an `nf_fields` lookup).
- `polarization_degrees` (Taelman's theorem; positive results need the existence hypotheses stated in the schema).
- The whole `k3_fq_lattices` computation (Oscar, Brandhorst–Muller).
- Positive realizability.  Only necessary conditions give −1 so far.

## Open decisions and next steps

1. **Cutoffs (waiting on David).**  All of F_2 is about 5.48M classes: 5,403,212 have ρ̄ ≥ 2 exactly, and the ρ̄ = 0 part is about 80k, not yet pinned down.  The sign −1 count is 1,485,591 (the notes said about 1.25M).  David asked about ρ̄ ≥ 16 for q = 3, ≥ 18 for q = 4, and ≥ 20 for q = 5, 7, 8.  That gives 7,868,568 rows, about 1.44 times F_2.  Changing q = 3 to ρ̄ ≥ 18 gives 5,599,956.  Cutoffs must be nondecreasing in q; the primitivity logic relies on this.  Exact counts for all q ≤ 64 are reproduced by `k3fq.py counts 3-64 --max-degree 6`.
2. **Lovelace run.**  `sh scripts/lovelace_counts.sh` gives the exact F_2 total.  All real data generation belongs on lovelace.  claudeF cannot SSH there, so David runs it.  mimir is shared and often heavily loaded (load 100–160 seen), so keep jobs there tiny and niced.
3. **Database test.**  Run `create_k3_fq_isog.py`, then `copy_from` the 12-row file, then `K3FqStats().setup()`.  Check real Postgres behaviour of `$startswith` on `poly_str`, `$overlaps` on `galois_groups`, and null ordering.  Add indexes on (q, poly) and on `poly_str` with `text_pattern_ops`.
4. **LMFDB integration.**
   - Create the `k3.fq.*` and `rcs.*.k3.fq` knowls from the drafts.  Several have [bracketed] placeholders for references, ranges and methods.
   - Add the k3fq link to the `motives.about` knowl (it renders `/Motive`).
   - Open the LMFDB PR when David says so.
5. **Mathematical conventions to confirm with David.**
   - "Igusa–Makai" in the notes was read as Gushel–Mukai.
   - The Verra fourfold point count assumes the 22-dimensional piece contains one of the three classes pulled back from P²×P².
   - Both are flagged in the schema notes.
6. **Performance at scale.**  The `basic` stage takes about 13 ms per polynomial.  The twists stage on full datasets is the likely bottleneck (it is pairwise within classes of at most 50).

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
