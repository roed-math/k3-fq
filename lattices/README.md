# Candidate Néron–Severi lattices with Frobenius (prototype)

For a class of K3 type over F_q with Weil polynomial

    P(t) = (1 - t)^ρ · C(t) · S(t),   C = ∏ Φ_d (d > 1),   E = lcm of the d,

these scripts list the pairs (N, φ) that could be (NS(X over F̄_q), Frobenius) for a K3 surface X in the class: N even of signature (1, ρ̄ − 1), φ an isometry with characteristic polynomial (x − 1)^ρ C(x), subject to the necessary conditions below.  This is a prototype for the `k3_fq_lattices` table, tested on q = 2, ρ̄ = 4.

## Method

Following Auel–Shaikh, *A cyclotomic Artin–Tate formula and sphere packing for K3 surfaces* (arXiv:2611.13468):

1. **Coinvariant part.** C_N = (N^φ)^⊥ with the restriction of φ is negative definite, even and has no roots (their Lemma 3.1), with characteristic polynomial C and order E.  By their Theorem 4.1 (m = E), |disc C_N| divides gcd(E^ρ, C(1)) · Res((t^E − 1)/(t − 1), S), and it is at least δ_r, the least discriminant of an even rootless lattice of rank r = ρ̄ − ρ (their Section 3).  For each admissible discriminant, Oscar's `enumerate_classes_of_lattices_with_isometry` lists the pairs (C_N, φ) genus by genus, and those with roots are discarded.
2. **Invariant part.** I = N^φ = NS(X) is even of signature (1, ρ − 1) with |disc I| = AT_1 / b² (Artin–Tate: |Br X| is a square).  Oscar's `representatives` lists each genus.
3. **Gluing.** Oscar's `equivariant_primitive_extensions` gives the pairs (N, φ) ⊃ (I, id) ⊕ (C_N, φ) up to the actions of O(I) and O(C_N, φ), so each isomorphism class of (N, φ) appears once.  The glue group has order dividing gcd(E^ρ, C(1)), which with Artin–Tate over F_{q^E} prunes the (I, C_N) pairs before gluing.
4. **Conditions on (N, φ)**, recorded as reasons for rejection:
   - `artin_tate`: for every m | E, AT_m / |disc N^{φ^m}| is a square integer (|Br X_m|).
   - `brauer_maps`: for m | m′ | E, the kernel of Br X_m → Br X_m′ is H¹(⟨φ^m⟩, N^{φ^{m′}}) (their Lemma 2.1), whose order is computed by their Lemma 1.3.  So it divides |Br X_m|, and |Br X_m| / |H¹| divides |Br X_m′|.
   - `ell_length`: for primes ℓ ≠ p, the ℓ-part of D(N) has at most 22 − ρ̄ cyclic factors.
   - `ell_glue`: for ℓ ≠ p, NS ⊗ Z_ℓ is primitive in the unimodular H²(X̄, Z_ℓ(1)), so D(N)_ℓ is glued Frobenius-equivariantly to the discriminant group of the transcendental lattice, which S(Frobenius) kills.  So S(φ) kills D(N) away from p.

   The coinvariant lattices of all the φ^m lie in C_N, so they are rootless too.  Nothing is checked at p yet (the crystalline condition).  `brauer_maps` and `ell_glue` go beyond the paper and should be checked by an expert.

Classes with E = 1 skip steps 1 and 3: N = I.

## Running

Julia 1.12 and Oscar 1.8 (versions pinned in `Manifest.toml`).  From the repository root:

    julia --project=lattices -e 'using Pkg; Pkg.instantiate()'
    sage -python lattices/prep.py TRANSFILE Q RHOBAR NSAMPLE SEED classes.jsonl
    julia --project=lattices lattices/k3lat.jl classes.jsonl lattices.jsonl

`prep.py` takes the transcendental parts of degree 22 − ρ̄ (one coefficient list per line, as written by `weilpoly`), a random sample of NSAMPLE of them (0 for all), times every cyclotomic part of degree ρ̄ containing 1 − t.  Over F_2 it adds the classes of the Kedlaya–Sutherland quartic census with that ρ̄.  It writes one JSON line per class: ρ, the orders d, E, the Artin–Tate values AT_m for m | E, the resultants of Theorem 4.1, and S.

`k3lat.jl` writes one JSON line per class: counts (`coinvariant_discs`, `coinvariant_lattices`, `candidates`, rejections by each condition, `lattices`) and the surviving pairs (Gram matrix and Frobenius matrix of N, acting on row vectors, and the discriminants of N, N^φ and C_N).  It skips classes already in the output, so it can be restarted.  Julia compiles the Oscar code paths on first use, which takes about 5 minutes per process; caches of coinvariant lattices, invariant lattices and gluings are kept per process.

## Test: q = 2, ρ̄ = 4 (November 2026, mimir)

200 random classes (20 transcendental parts × the 10 cyclotomic parts with ρ ≥ 1) and 100 classes from the quartic census:

- **Census:** all 109 census classes in the run (the 100, plus 9 among the random ones) keep at least one pair (N, φ), as they must.
- **Exclusions:** of the 108 random classes that pass the basic necessary conditions of `k3_fq_isog`, 55 have no pair left.  32 already fail at step 1 (no rootless coinvariant lattice with an admissible discriminant).  The other 23 fail the conditions in step 4 after gluing.
- **E ≤ 2:** 52 of those 108 classes, 21 excluded, 41 pairs for the other 31.
- **Pairs:** at most 4 per class.  All 399 output pairs were checked independently in Sage (φ an isometry with the right characteristic polynomial, invariant rank ρ, signature (1, 3)).
- **Time:** after compilation, a median of 0.02 s per class and about 0.7 s on average once the caches fill, with a maximum of 71 s (classes with Φ_2Φ_3 and about 30 admissible coinvariant discriminants).  All 183,780 classes with ρ̄ = 4 over F_2 would take very roughly 25–35 CPU-hours.

## Known issues

- Oscar's `root_test = true` option of `enumerate_classes_of_lattices_with_isometry` skips the genus of A₂(−2) with Φ₃ (minimum −4, so rootless), so the rootless filter is done here instead.
- Possible speedups: apply Theorem 4.1 for every m | E to the coinvariant lattice before gluing (most rejected candidates fail `artin_tate` at an intermediate m), share caches across processes, and build a Julia system image to avoid recompiling.
