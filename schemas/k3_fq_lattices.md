Table name: `k3_fq_lattices`

This table stores the candidate geometric Néron–Severi lattices, with their Frobenius action, for the classes in `k3_fq_isog`.

Let $X/\mathbb{F}_q$ be a K3 surface whose class in `k3_fq_isog` has polynomial $P$, with cyclotomic part $C$ of degree $\bar\rho$.  Then $N = NS(X_{\overline{\mathbb{F}}_q})$ is an even lattice of signature $(1, \bar\rho - 1)$.  Frobenius acts on it by an isometry $\varphi$ of finite order with characteristic polynomial $\pm C(T)$ (reversed).  The action preserves the ample cone, which is a Weyl chamber of the positive cone.  The invariant sublattice $N^{\varphi}$ is $NS(X)$.

In some cases $P$ determines $(N, \varphi)$ up to isomorphism.  In general there are finitely many possibilities, and the rows of this table enumerate them.  The enumeration uses Oscar's lattices-with-isometry machinery (Brandhorst–Muller) and is subject to the necessary conditions listed under `conditions`.  So each row is an isomorphism class of pairs $(N, \varphi)$ that is not ruled out for $P$; it need not come from an actual K3 surface.

**Label**: `{isog_label}.{n}` where `n` is a counter (starting at 1) ordering the pairs by: `genus_label`, then `invariant_disc` (absolute value), then the canonical Gram matrix and Frobenius matrix.

| Column | Type | Description |
| --- | --- | --- |
| label | text | label `{isog_label}.{n}` |
| isog_label | text | label of the class in `k3_fq_isog` |
| q | integer | cardinality of the base field (duplicated from `k3_fq_isog` for searching) |
| p | integer | characteristic of the base field |
| rank | smallint | $\bar\rho$, the rank of $N$ (the geometric Picard rank of the class) |
| gram | integer[] | Gram matrix of $N$ (row by row, $\bar\rho^2$ entries), in a reduced basis |
| frobenius | integer[] | matrix of $\varphi$ in the basis of `gram`, row by row, acting on row vectors: $v \mapsto vF$.  It satisfies $F\,G\,F^{T} = G$ |
| frobenius_order | integer | order of $\varphi$; this equals `geometric_extension_degree` of the class |
| genus_label | text | label of the genus of $N$ in `lat_genera` |
| lattice_label | text | label of $N$ in `lat_lattices`, if present |
| disc | bigint | determinant of `gram` (sign $(-1)^{\bar\rho - 1}$) |
| discriminant_group_invs | integer[] | invariants of the discriminant group $N^\vee/N$ |
| invariant_rank | smallint | rank of $N^{\varphi}$, equal to `picard_rank` of the class |
| invariant_gram | integer[] | Gram matrix of $N^{\varphi} = NS(X)$ in a reduced basis |
| invariant_disc | bigint | determinant of `invariant_gram` |
| invariant_genus_label | text | label of the genus of $N^{\varphi}$ in `lat_genera` |
| brauer_order | numeric | $AT_1/\lvert\mathrm{disc}\, N^{\varphi}\rvert$, where $AT_1$ = `k3_fq_isog.artin_tate_value`.  By the Artin–Tate formula this is $\#\mathrm{Br}(X)$; it must be a perfect square |
| brauer_orders | numeric[] | $AT_n/\lvert\mathrm{disc}\, N^{\varphi^n}\rvert$ for $n$ in `k3_fq_isog.extension_degrees`: the orders of $\mathrm{Br}(X_{\mathbb{F}_{q^n}})$ |
| min_quasi_ample_degree | integer | $\min\{v^2 : v \in N^\varphi,\ v^2 > 0\}$: the smallest degree of a big and nef class defined over $\mathbb{F}_q$.  Any positive class can be moved into the nef cone by reflections, so this is the minimal positive norm in $N^\varphi$ |
| min_ample_degree | integer | the smallest $v^2$ for $v \in N^{\varphi}$ in the interior of a $\varphi$-stable Weyl chamber (so $v^2 > 0$ and $v \cdot r \ne 0$ for every $r \in N$ with $r^2 = -2$), minimized over $\varphi$-stable chambers.  Null if $\bar\rho = 0$ |
| has_roots | boolean | whether $N$ contains a vector $r$ with $r^2 = -2$, i.e. whether $X_{\overline{\mathbb{F}}_q}$ contains smooth rational curves |
| contains_kummer_lattice | boolean | whether $N$ contains the Kummer lattice (rank 16, discriminant $2^6$) primitively.  By Nikulin's criterion this is the lattice-theoretic condition for $X_{\overline{\mathbb{F}}_q}$ to be a Kummer surface |
| conditions | text[] | the necessary conditions checked for this pair.  Possible entries: `artin_tate` (`brauer_orders` are squares), `ell_adic` (for each $\ell \ne p$, $N \otimes \mathbb{Z}_\ell$ embeds primitively in the even unimodular $\mathbb{Z}_\ell$-lattice of rank 22, with complement carrying the transcendental part), `crystalline` (compatibility with the $p$-adic Newton polygon: Ogus' conditions on $N \otimes \mathbb{Z}_p$), `chamber` ($\varphi$ preserves a Weyl chamber) |
| realizable_k3 | smallint | 1 if there is a K3 surface over $\mathbb{F}_q$ in this class with $(NS(X_{\overline{\mathbb{F}}_q}), \mathrm{Frob}) \cong (N, \varphi)$, $-1$ if not, 0 if unknown |
| realizable_nck3 | smallint | as `realizable_k3`, for noncommutative K3 surfaces |
| realizable_kummer | smallint | as `realizable_k3`, for Kummer surfaces |

### Notes

- For ordinary classes one can also work with the full lattice: by Taelman, an ordinary K3 surface is determined by $H^2$ of its canonical lift, an even unimodular lattice $\Lambda$ of signature $(3,19)$, with a Frobenius $F$ and a choice of ample cone.  Here $F$ is a similitude, $\langle Fx, Fy\rangle = q^2\langle x, y\rangle$, so $F/q$ is a $\mathbb{Z}[1/p]$-isometry whose reciprocal characteristic polynomial is $P(T)/q$.  (The existence half of Taelman's theorem assumes potential semistable reduction; see `k3_fq_isog.realizable_k3`.)  If we decide to store these *Deligne–Taelman modules*, they should go in a third table keyed by `isog_label`, with the Frobenius matrix on $\mathbb{Z}^{22}$ in place of `frobenius`.
- The enumeration for pointless classes over $\mathbb{F}_2$ (Brandhorst) produced about a million genera and took about 100 CPU years, so `num_lattices` in `k3_fq_isog` will be null for most classes for some time.
- This table only concerns the K3 interpretation.  Lattice data for the fourfold interpretations is different: for a cubic fourfold it is the positive definite lattice of algebraic classes in the full $H^4(Y,\mathbb{Z})$, containing $h^2$ of norm $3$, together with its Hassett markings.  Relating it to a K3 or Mukai lattice needs an explicit construction (e.g. Hassett's associated K3 surfaces for admissible discriminants).  If wanted, that data belongs in separate tables, with their own realizability columns.
