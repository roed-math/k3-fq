Table name: `k3_fq_isog`

This table stores motives of K3 type over finite fields up to isogeny: equivalently, their zeta functions.  Each row is a polynomial
$$P(T) = a_0 + a_1 T + \cdots + a_{22} T^{22} \in \mathbb{Z}[T]$$
with $a_0 = q$, $a_{22} = \varepsilon q$ for a sign $\varepsilon \in \{\pm 1\}$, and all complex roots on the unit circle.  Here $q = p^r$ is a prime power.  The coefficients satisfy $a_{22-i} = \varepsilon a_i$.

If $X/\mathbb{F}_q$ is a K3 surface then $P(T) = q \det(1 - q^{-1}\mathrm{Frob}\, T \mid H^2_{et}(X_{\overline{\mathbb{F}}_q}, \mathbb{Q}_\ell))$, and $\varepsilon$ is the sign of the functional equation.  We write $P(T) = q\prod_{i=1}^{22}(1 - \alpha_i T)$; the $\alpha_i$ are the *normalized Frobenius eigenvalues* (Frobenius eigenvalues divided by $q$).  They lie on the unit circle and $\prod \alpha_i = \varepsilon$.  The usual L-polynomial is $L(T) = \det(1 - \mathrm{Frob}\, T \mid H^2) = q^{-1}P(qT) = 1 + \sum_{i \ge 1} q^{i-1}a_i T^i$.  Every $L$ for a K3 surface arises this way: by Mazur's theorem, $q^{i-1} \mid$ the $i$th coefficient of $L$.

The *cyclotomic part* $C(T)$ of $P$ is the product of its cyclotomic factors (normalized with $C(0) = 1$).  The *transcendental part* is $P_T = P/C \in \mathbb{Z}[T]$, which has constant term $q$ and no roots of unity among its roots.  By the Tate conjecture (a theorem for K3 surfaces), the roots of $C$ correspond to algebraic classes.  So $\deg C$ is the geometric Picard rank.  The multiplicity of $T = 1$ as a root of $P$ is the Picard rank over $\mathbb{F}_q$.

Other varieties of K3 type are treated through the same 22-dimensional motive $M$ of weight $0$, with $P(T) = q\det(1 - \mathrm{Frob}\,T \mid M)$.  For each type, $M$ is the following piece of the cohomology:
- K3 surface: $H^2(X)(1)$.
- cubic fourfold: $H^4_{prim}(Y)(2)$.
- Gushel–Mukai fourfold: $H^4_{van}(Y)(2)$.
- Verra fourfold: the K3-type piece of $H^4(Y)(2)$ (see the notes below).

**Label**: `q.ρ.c.t` where

- `q` is the cardinality of the base field;
- `ρ` is the geometric Picard rank $\bar\rho$, which is even, and $k = (22 - \bar\rho)/2$;
- `c` encodes the cyclotomic part $C = (1-T)^{e_1}\prod_{m > 1}\Phi_m^{e_m}$.  Index the exponent vectors $(e_1, e_2, e_3, \ldots)$ by the 43 integers $m$ with $\varphi(m) \le 22$, in increasing order.  Sort the vectors with $\sum_m e_m\varphi(m) = \bar\rho$ in decreasing lexicographic order and number them from $0$.  Then `c` is the number of $C$, written in base 26 with digits `a`–`z` (`a` = 0).  So `a` is $(1-T)^{\bar\rho}$, and for $\bar\rho = 2$ the codes `a`–`f` are $(1-T)^2$, $(1-T)(1+T)$, $(1+T)^2$, $\Phi_3$, $\Phi_4$, $\Phi_6$.  There are 60,280 cyclotomic parts of degree 22, so `c` has at most four letters;
- `t` encodes the coefficients $b_1, \ldots, b_k$ of the transcendental part $P_T = q + b_1 T + \cdots + b_{2k-1}T^{2k-1} + qT^{2k}$.  Each $b_i$ is written in balanced base 26, where the letters `a`–`z` stand for the digits $-13, \ldots, 12$ (so `n` = 0, `o` = 1, `m` = $-1$, and `oa` = 13); no sign is needed.  The codes are padded with leading `n`s to a common width $w$, the smallest that fits all of them, and concatenated, so `t` has $kw$ letters.  Usually $w = 1$: over $\mathbb{F}_2$ this holds for more than 98% of the transcendental parts of each degree up to 16.

The transcendental part is palindromic ($b_{2k-i} = b_i$), so the label determines $P = C \cdot P_T$.  When $\bar\rho = 22$ there is no transcendental part, and the label is `q.22.c`.  The sign $\varepsilon = (-1)^{e_1}$ is determined by `c`.  For example, $2(1-T)^{22}$ has label `2.22.a` and $2(1-T)^{21}(1+T)$ has label `2.22.b`.  Over $\mathbb{F}_4$, $(1-T)^{16}\,\Phi_4\,(4 + T + 6T^2 + T^3 + 4T^4)$ has label `4.18.e.ot`, and $(1-T)^{18}(2 - 3T + 2T^2)^2$ has label `4.18.a.nboe` ($b_1 = -12$, $b_2 = 17$).

Labels are ordered by $q$, then $\bar\rho$, then the number encoded by `c`, then $b_1, \ldots, b_k$ as integers.  This order defines the representative in `geometric_class`.  They are implemented in `k3fq/labels.py` in [roed-math/k3-fq](https://github.com/roed-math/k3-fq).

| Column | Type | Description |
| --- | --- | --- |
| label | text | label `q.ρ.c.t`, as above |
| q | integer | cardinality of the base field |
| p | integer | characteristic of the base field |
| r | smallint | degree of the base field over its prime field, so that $q = p^r$ |
| sign | smallint | sign $\varepsilon \in \{1, -1\}$ of the functional equation; equivalently $(-1)^{\rho}$ where $\rho$ is the multiplicity of $T=1$ as a root of $P$.  When $\varepsilon = -1$, $(1-T)(1+T)$ divides $P$ |
| poly | bigint[] | coefficients $[a_0, a_1, \ldots, a_{22}]$ of $P(T)$, with $a_0 = q$ and $a_{22} = \varepsilon q$ |
| poly_str | text | space-separated string of the coefficients of `poly`, with a trailing space, for searching on initial segments |
| real_poly | bigint[] | coefficients (constant term first) of the polynomial $h$ with leading coefficient $q$ and roots $\alpha_j + \bar\alpha_j = 2\cos(\arg \alpha_j)$.  For $\varepsilon = 1$ it is defined by $P(T) = T^{11} h(T + T^{-1})$ and has degree 11.  For $\varepsilon = -1$ it is defined by $P(T) = (1-T^2)\,T^{10} h(T + T^{-1})$ and has degree 10 |
| angles | double precision[] | arguments of the roots of $P$ in the closed upper half plane, divided by $\pi$, sorted, with multiplicity (as for abelian varieties).  All lie in $[0,1]$; a root $\pm 1$ of multiplicity $m$ contributes $m$ copies of $0$ or $1$ |
| angle_rank | smallint | the dimension of the $\mathbb{Q}$-span of $\{\arg(\alpha_i)/2\pi\} \cup \{1\}$, minus 1.  This equals the rank of the multiplicative group generated by the $\alpha_i$, and also the dimension of the Frobenius torus.  It lies between $0$ and $(22 - \bar\rho)/2$, and is $0$ exactly when $\bar\rho = 22$ |
| angle_corank | smallint | $(22 - \bar\rho)/2$ minus `angle_rank`, where $\bar\rho$ = `geometric_picard_rank`: the number of multiplicative relations among the transcendental roots |
| slopes | text[] | the $q$-adic Newton slopes of $L(T)$ (equivalently of Frobenius on $H^2$), as strings of rational numbers in $[0, 2]$.  The list is sorted, with repeated slopes suffixed `A`, `B`, ... as in `av_fq_isog` |
| height | smallint | the (Artin–Mazur) height $h$: the number of slopes less than $1$, counted with multiplicity.  For a K3 surface this is the height of the formal Brauer group, and the slopes are $1 - 1/h$ ($h$ times), $1$ ($22 - 2h$ times) and $1 + 1/h$ ($h$ times).  It is null when supersingular (height $\infty$), and $1$ when ordinary.  When $r > 1$ there are polynomials whose Newton polygon is not of this form (e.g. $P = 2P'$ with $P'$ of K3 type over $\mathbb{F}_2$, viewed over $\mathbb{F}_4$); for these `has_k3_newton_polygon` is false |
| has_k3_newton_polygon | boolean | whether the slopes have the shape $1 - 1/h, 1, \ldots, 1, 1 + 1/h$ for some $1 \le h \le 11$, or are all equal to $1$.  This is automatic when $q$ is prime |
| newton_elevation | smallint | analogous to `av_fq_isog.newton_elevation`.  It is the number of lattice points $(x,y)$ with $0 \le x \le 11$ that lie on or above the ordinary Newton polygon (vertices $(0,0),(1,0),(21,20),(22,22)$) and strictly below the Newton polygon of $L$.  For K3-type Newton polygons it equals $h-1$, and it is $11$ in the supersingular case.  Null when `has_k3_newton_polygon` is false |
| newton_coelevation | smallint | the number of lattice points $(x,y)$ with $0 \le x \le 11$ lying on or above the Newton polygon of $L$ and strictly below the supersingular polygon $y = x$.  For K3-type Newton polygons it equals $12 - h$ and is $0$ when supersingular.  So `newton_elevation + newton_coelevation = 11`.  Null when `has_k3_newton_polygon` is false |
| is_ordinary | boolean | whether $p \nmid a_1$, i.e. slope $0$ occurs.  Then the Newton polygon is of K3 type with height $1$.  (Over non-prime fields, `height` $= 1$ alone does not imply this: e.g. $2(2 - T + 2T^2)(1-T)^{20}$ over $\mathbb{F}_4$ has slopes $1/2, 1, \ldots, 1, 3/2$) |
| is_supersingular | boolean | whether all slopes equal $1$.  Equivalently all roots of $P$ are roots of unity, i.e. `geometric_picard_rank` = 22 |
| cyclotomic_factors | integer[] | the sorted list, with multiplicity, of the $m$ such that $C(T) = \pm\prod \Phi_m(T)$ |
| picard_rank | smallint | $\rho$, the multiplicity of $T = 1$ as a root of $P$.  For a K3 surface this is the rank of $NS(X)$ |
| geometric_picard_rank | smallint | $\bar\rho = \deg C$, the number of roots of $P$ that are roots of unity.  For a K3 surface this is the rank of $NS(X_{\overline{\mathbb{F}}_q})$ |
| geometric_extension_degree | integer | the least common multiple $E$ of the entries of `cyclotomic_factors` (1 if there are none).  It is the smallest $n$ with $\rho(X_{\mathbb{F}_{q^n}}) = \bar\rho$ |
| extension_degrees | integer[] | the sorted list $D$ of least common multiples of subsets of `cyclotomic_factors` (including $1$, the empty lcm).  These are the degrees $n$ at which the Picard rank jumps: for any $n$, $\rho_n = \rho_{n'}$ where $n' = \mathrm{lcm}\{m : m \mid n\} \in D$ |
| picard_ranks | smallint[] | the Picard ranks $\rho_n = \sum_{m \mid n} \#\{m \text{ in cyclotomic\_factors}\}\cdot\varphi(m)$ for $n$ in `extension_degrees` |
| artin_tate_values | numeric[] | for $n$ in `extension_degrees`, the positive integer $AT_n = \left(P_n(T)/(1-T)^{\rho_n}\right)\big|_{T=1}$.  Here $P_n(T) = q^n\prod(1 - \alpha_i^n T)$ is the corresponding polynomial over $\mathbb{F}_{q^n}$.  For a K3 surface, the Artin–Tate formula gives $AT_n = \#\mathrm{Br}(X_{\mathbb{F}_{q^n}}) \cdot \lvert\mathrm{disc}\, NS(X_{\mathbb{F}_{q^n}})\rvert$ |
| artin_tate_value | numeric | $AT_1$, duplicated for searching |
| ns_disc_class | bigint | the squarefree integer $(-1)^{\rho-1}\,\mathrm{sqf}(AT_1)$.  For a K3 surface this is the class of $\mathrm{disc}\,NS(X)$ in $\mathbb{Q}^\times/\mathbb{Q}^{\times 2}$, since the order of the Brauer group is a square.  Null if $\rho = 0$ |
| ej_condition | boolean | whether $P$ satisfies the Artin–Tate condition of Elsenhans and Jahnel, as tested in Kedlaya–Sutherland (`k3-scripts/ej-test.sage` in root-unitary).  Write $S = P/(1-T)^{\rho}$; the condition is that $S(-1)$ is a perfect square ($0$ allowed).  When $(1+T) \nmid P$ this says that $AT_1 AT_2$ is a square, which holds for K3 surfaces |
| transcendental_poly | bigint[] | coefficients of $P_T = P/C$, constant term $q$; just `[q]` when supersingular |
| transcendental_factors | jsonb | factorization of $P_T$ over $\mathbb{Q}$: a list of pairs `[f, e]`, with each `f` the coefficient list (constant term first) of a primitive irreducible factor with positive constant term, and `e` its multiplicity.  When $q$ is prime this always has a single factor with $e = 1$ |
| transcendental_irreducible | boolean | whether $P_T$ is irreducible over $\mathbb{Q}$ (as opposed to a power of an irreducible polynomial, or a product of distinct ones; this needs $r > 1$).  Null when supersingular |
| number_fields | text[] | for each factor in `transcendental_factors`, the LMFDB label of the CM field $\mathbb{Q}[T]/(f)$, or the empty string if it is not in the LMFDB |
| number_field_coeffs | jsonb | for each factor in `transcendental_factors`, the coefficients of a reduced defining polynomial for $\mathbb{Q}[T]/(f)$.  This is `polredabs` when it finishes within the time limit, otherwise `polredbest` |
| galois_groups | text[] | for each factor in `transcendental_factors`, the label `nTk` of the Galois group of $f$ |
| k3_counts | numeric[] | $\#X(\mathbb{F}_{q^n}) = 1 + q^{2n} + q^n\sum_i \alpha_i^n$ for $n = 1, \ldots, 10$: the point counts of a K3 surface in this class |
| k3_count | bigint | $\#X(\mathbb{F}_q)$, duplicated for searching (equal to $1 + q^2 - a_1$) |
| cubic_counts | numeric[] | $1 + q^n + q^{2n} + q^{3n} + q^{4n} + q^{2n}\sum_i\alpha_i^n$ for $n = 1,\ldots,10$: the point counts of a cubic fourfold $Y$ with $H^4_{prim}(Y)(2)$ in this class |
| cubic_count | bigint | the first entry of `cubic_counts` (the number of $\mathbb{F}_q$-points of a cubic fourfold in this class), duplicated for searching |
| gm_counts | numeric[] | $1 + q^n + 2q^{2n} + q^{3n} + q^{4n} + q^{2n}\sum_i\alpha_i^n$ for $n = 1,\ldots,10$: the point counts of a Gushel–Mukai fourfold with vanishing cohomology in this class |
| gm_count | bigint | the first entry of `gm_counts` (the number of $\mathbb{F}_q$-points of a Gushel–Mukai fourfold in this class), duplicated for searching |
| verra_counts | numeric[] | $1 + 2q^n + 2q^{2n} + 2q^{3n} + q^{4n} + q^{2n}\sum_i\alpha_i^n$ for $n = 1,\ldots,10$: the point counts of a Verra fourfold whose K3-type piece is in this class (see the notes below) |
| verra_count | bigint | the first entry of `verra_counts` (the number of $\mathbb{F}_q$-points of a Verra fourfold in this class), duplicated for searching |
| has_nonnegative_counts | boolean | whether $\#X(\mathbb{F}_{q^n}) \ge 0$ for all $n \ge 1$ (as a K3 surface).  This only needs checking for the finitely many $n$ with $q^n < 22$ |
| has_count_growth | boolean | whether the number of closed points of degree $n$, $\frac{1}{n}\sum_{d \mid n}\mu(n/d)\#X(\mathbb{F}_{q^d})$, is nonnegative for every $n \ge 1$ (as a K3 surface).  This implies $\#X(\mathbb{F}_{q^m}) \ge \#X(\mathbb{F}_{q^n})$ whenever $n \mid m$, e.g. more points over $\mathbb{F}_{16}$ than over $\mathbb{F}_4$ |
| is_pointless | boolean | whether `k3_count` is $0$ |
| is_primitive | boolean | whether this class is not the base change of a class over a proper subfield of $\mathbb{F}_q$ |
| primitive_models | text[] | labels of classes over proper subfields $\mathbb{F}_{q_0}$ ($q_0^d = q$) whose base change to $\mathbb{F}_q$ is this class.  Null if primitive |
| geometric_class | text | label of the representative of this class's twist class.  Two classes over $\mathbb{F}_q$ are *twists* if they become equal over some finite extension, i.e. $P_n = P'_n$ for some $n$.  The representative is the member of the twist class with the smallest label, in the order defined above (all members have the same $q$ and $\bar\rho$).  The class of $(1-T)^{\bar\rho}P_T$ is always a twist, so when all classes with this $q$ and $\bar\rho$ are present (as in the database), the representative has cyclotomic code `a` |
| twist_count | integer | number of classes in the twist class (including this one) |
| twists | jsonb | list of triples `[label, bc_label, n]`, where `label` is a twist, `n` is a minimal degree (for divisibility) over which the two classes become equal, and `bc_label` is the label of the common base change to $\mathbb{F}_{q^n}$ (null when $n$ times the number of bits of $q$ exceeds $128$).  This is the same convention as `av_fq_isog.twists`.  Null when `twist_count` exceeds 50 (supersingular classes and classes with large cyclotomic parts have very large twist classes); search on `geometric_class` instead |
| max_twist_degree | integer | the largest `n` appearing in `twists` (null when `twists` is) |
| kummer_av_classes | text[] | labels of isogeny classes of abelian surfaces $A/\mathbb{F}_q$ (in `av_fq_isog`) whose wedge square, normalized, has the same transcendental part as $P$.  A Kummer surface $\mathrm{Km}(A)$ ($p > 2$) has $P = (\wedge^2$ part$) \cdot (\text{permutation part on } A[2])$, so these are the candidates for Kummer structures.  Empty if none, null if not computed |
| realizable_k3 | smallint | 1 if this class is known to contain a K3 surface over $\mathbb{F}_q$, $-1$ if it is known not to, 0 otherwise.  Necessary conditions include $\rho \ge 1$, `ej_condition` and `has_count_growth`.  A positive value must come with a proof of existence; at present these are the classes over $\mathbb{F}_2$ that contain a smooth quartic surface, from the census of Kedlaya and Sutherland.  For ordinary classes, Taelman's description of ordinary K3 surfaces is fully faithful unconditionally, but essential surjectivity assumes potential semistable reduction.  That is known in some cases (e.g. Taelman's Theorem 4.6: geometric Picard rank $\ge 12$ and $p \ge 5$), so existence in characteristics $2$ and $3$ needs separate justification |
| realizable_nck3 | smallint | same for noncommutative K3 surfaces (K3 categories over $\mathbb{F}_q$, e.g. twisted K3 surfaces and Kuznetsov components), with $M$ the K3-type part of the Mukai motive |
| realizable_cubic | smallint | same for cubic fourfolds, with $M = H^4_{prim}(Y)(2)$.  The relevant lattice is the positive definite lattice $A(Y)$ of algebraic classes in the full $H^4(Y, \mathbb{Z})$, which contains $h^2$ (norm $3$).  Hassett markings are saturated rank 2 sublattices $\langle h^2, v\rangle$ with Gram matrix $\begin{pmatrix}3 & a\\ a & b\end{pmatrix}$ and discriminant $d = 3b - a^2$.  The algebraic classes in $M$ form the even lattice $A(Y) \cap (h^2)^\perp$, which is not the Néron–Severi lattice of a K3 surface (see `k3_fq_lattices`) |
| realizable_nccubic | smallint | same for Kuznetsov components of cubic fourfolds |
| realizable_gm | smallint | same for Gushel–Mukai fourfolds |
| realizable_ncgm | smallint | same for Kuznetsov components of Gushel–Mukai fourfolds |
| realizable_verra | smallint | same for Verra fourfolds |
| realizable_ncverra | smallint | same for Kuznetsov components of Verra fourfolds |
| realizable_kummer | smallint | same for Kummer surfaces $\mathrm{Km}(A)$ over $\mathbb{F}_q$.  It is $-1$ if `kummer_av_classes` is empty or `realizable_k3` is $-1$.  Null in characteristic 2 |
| polarization_degrees | integer[] | the sorted degrees $2d \le 100$ for which the class is known to contain a K3 surface over $\mathbb{F}_q$ with a primitive polarization of degree $2d$ defined over $\mathbb{F}_q$.  At present they come from the quartic census over $\mathbb{F}_2$ (degree $4$); Taelman's theorem would give them for ordinary classes.  Degrees whose existence depends on potential semistable reduction (such as Taelman's results for $p \le 3$) are not included.  Null if not computed |
| polarization_degrees_complete | boolean | whether `polarization_degrees` is known to contain every degree $2d \le 100$ that occurs, so that the other degrees do not.  False for the lists from the census (for instance, a K3 surface with a primitive polarization of degree $4$ need not be a smooth quartic).  Null when `polarization_degrees` is null |
| num_lattices | integer | number of rows in `k3_fq_lattices` for this class (null if the lattice computation has not been run for this class) |
| min_ample_degree | integer | minimum of `k3_fq_lattices.min_ample_degree` over the lattices for this class |
| min_quasi_ample_degree | integer | minimum of `k3_fq_lattices.min_quasi_ample_degree` over the lattices for this class |

### Notes

- *Sign $-1$.*  The enumeration for $\varepsilon = -1$ uses $P = (1 - T^2)P'$ with $P'$ of degree 20 and sign $+1$, exactly as `WeilPolynomials(..., sign=-1)` does internally.  Note that $\varepsilon = -1$ forces $\rho$ to be odd, and $\varepsilon = 1$ forces it to be even.
- *Point counts for other types.*  The formulas above assume that the cohomology outside $M$ consists of Tate classes defined over $\mathbb{F}_q$.  They have Betti numbers $b_2 = 1, b_4 = 23$ for cubic fourfolds and $b_2 = 1, b_4 = 24$ for Gushel–Mukai fourfolds.  A Verra fourfold (double cover of $\mathbb{P}^2\times\mathbb{P}^2$ branched along a $(2,2)$ divisor) has $\chi = 30$, $b_2 = b_6 = 2$, $b_4 = 24$, with three algebraic classes pulled back from $\mathbb{P}^2\times\mathbb{P}^2$ in $H^4$.  The formula above takes the 22-dimensional piece to contain one of them, leaving two Tate classes outside.  This should be checked against the intended definition of "Verra fourfold of K3 type".
- *"Igusa–Makai"* in the discussion notes has been read as *Gushel–Mukai*.
- *Twists.*  The geometric class is computed from the multiset of $\{M\theta_j\} \bmod 1$ over the transcendental angles $\theta_j = \arg(\alpha_j)/2\pi$, where $M = \mathrm{lcm}\{m : \varphi(m) \le 44\}$.  The minimal degrees are then verified by exact base change.
- *Not yet computed* (left null by the current pipeline): `number_fields` labels (pending a lookup against `nf_fields`), `realizable_*` beyond the necessary conditions and the quartic census, `polarization_degrees` beyond the census, and the lattice summary columns.
- *Possible further columns*: the rational endomorphism algebra of the motive (center $\mathbb{Q}[T]/(f)$ for the transcendental part, with local invariants as in `av_fq_endalg_data`), and the Mumford–Tate group (the Frobenius torus has dimension `angle_rank`).
