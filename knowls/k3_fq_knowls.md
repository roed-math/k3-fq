# Draft knowls for motives of K3 type over finite fields (`k3.fq.*`)

These are drafts for the knowls referenced by the `k3fq` section of the LMFDB
(`lmfdb/k3fq/` on the `k3fq` branch, url `/Motive/K3/Fq/`).  Each section gives
the knowl id, a proposed title, and the proposed content (LMFDB knowl markdown,
with `{{ KNOWL(...) }}` cross references and `\cite{...}` as in existing knowls).
Text in [square brackets] marks places that need input.

Existing knowls that the section reuses (unchanged): `ag.base_field`,
`intro.search`, `nf.galois_group.name`, and the column/statistics knowls
`stats.buckets`, `stats.totals`, `stats.proportions` used by every dynamic
statistics page.

Knowls for interesting classes follow the usual pattern `k3.fq.<label>.top`
(category `k3`), and are listed on `/Motive/K3/Fq/interesting`.

Suggested edit to an existing knowl: the page `/Motive` renders the knowl
`motives.about`, which currently lists only hypergeometric motives.  Add

    <li><a href="{{ url_for('k3fq.index') }}">Motives of K3 type over finite fields</a>

after the hypergeometric line.

---

## k3.fq.lmfdb_label

**Title:** Labels for motives of K3 type over finite fields

The label of a class of {{ KNOWL('k3.fq.k3_polynomial', 'motives of K3 type') }} over $\F_q$ has the form **q.s.c1_c2_..._c11**, where

- **q** is the cardinality of the base field $\F_q$;
- **s** is `p` if the {{ KNOWL('k3.fq.sign', 'sign') }} $\varepsilon$ is $+1$, and `m` if $\varepsilon = -1$;
- **c1, ..., c11** encode the coefficients $a_1, \ldots, a_{11}$ of the polynomial
$$P(T) = q + a_1 T + a_2 T^2 + \cdots + a_{22} T^{22},$$
as in the {{ KNOWL('av.fq.lmfdb_label', 'labels of isogeny classes of abelian varieties') }}: nonnegative integers are written in base 26 with the digits a, b, ..., z (where a = 0), and negative integers get a leading a.

Since $a_{22-i} = \varepsilon a_i$, the label determines $P(T)$.  When $\varepsilon = -1$ the middle coefficient $a_{11}$ vanishes, so the last code is always `a`.  The sign is part of the label since a polynomial with sign $+1$ can also have $a_{11} = 0$.

For example, $P(T) = 2(1-T)^{22}$, which is supersingular over $\F_2$, has $a_1 = -44$, encoded as `abs`, and label `2.p.abs_ru_aeom_vqs_aczxs_imtm_atkou_bkkbs_aceprc_cvpea_adchca`.

The {{ KNOWL('k3.fq.lattices', 'candidate Néron–Severi lattices') }} of a class with label $L$ have labels $L.n$, where $n = 1, 2, \ldots$ is a counter.

## k3.fq.k3_polynomial

**Title:** Motive of K3 type over a finite field

A **motive of K3 type** over a finite field $\F_q$ is a 22-dimensional motive $M$ of weight $0$ that looks like the second cohomology of a K3 surface: for a K3 surface $X/\F_q$ it is $H^2(X)(1)$; for a cubic fourfold $Y$ it is $H^4_{\mathrm{prim}}(Y)(1)$; for a Gushel–Mukai fourfold it is $H^4_{\mathrm{van}}(Y)(1)$; and for a Verra fourfold it is the K3-type piece of $H^4(Y)(1)$.

The LMFDB classifies such motives by the polynomial
$$P(T) = q \det\left(1 - q^{-1}\mathrm{Frob}\, T \mid H^2_{\text{ét}}(X_{\overline{\F}_q}, \Q_\ell)\right) = q\prod_{i=1}^{22}(1 - \alpha_i T) = q + a_1 T + \cdots + a_{22}T^{22} \in \Z[T],$$
where the **normalized Frobenius eigenvalues** $\alpha_i$ (the eigenvalues of Frobenius divided by $q$) lie on the unit circle.  So $a_0 = q$, $a_{22} = \varepsilon q$ and $a_{22-i} = \varepsilon a_i$, where $\varepsilon = \prod \alpha_i = \pm 1$ is the {{ KNOWL('k3.fq.sign', 'sign') }}.

The database contains every polynomial of this shape for the base fields listed in the {{ KNOWL('rcs.cande.k3.fq', 'completeness') }} knowl, whether or not it comes from an actual K3 surface or fourfold; see {{ KNOWL('k3.fq.realizability', 'realizability') }}.  Each such polynomial is called a *class*.  The usual {{ KNOWL('k3.fq.l_polynomial', 'L-polynomial') }} is $L(T) = q^{-1}P(qT)$.

## k3.fq.l_polynomial

**Title:** L-polynomial of a motive of K3 type

The **L-polynomial** of a {{ KNOWL('k3.fq.k3_polynomial', 'motive of K3 type') }} over $\F_q$ is the characteristic polynomial of Frobenius on $H^2$,
$$L(T) = \det(1 - \mathrm{Frob}\, T \mid H^2) = q^{-1} P(qT) = 1 + \sum_{i=1}^{22} q^{i-1} a_i T^i,$$
where $P(T) = q + a_1 T + \cdots + a_{22} T^{22}$ is the normalized polynomial used in the LMFDB.  For a K3 surface $X/\F_q$ the zeta function is
$$Z(X, T) = \frac{1}{(1-T)\, L(T)\, (1 - q^2 T)}.$$

Conversely, every polynomial $L(T)$ with constant term $1$ whose reciprocal roots have absolute value $q$ and with $q^{i-1}$ dividing the coefficient of $T^i$ arises from a polynomial $P(T)$ in this way; the divisibility holds for K3 surfaces by Mazur's theorem.

## k3.fq.sign

**Title:** Sign of a motive of K3 type

The **sign** of a {{ KNOWL('k3.fq.k3_polynomial', 'motive of K3 type') }} with polynomial $P(T) = q\prod_{i=1}^{22}(1-\alpha_i T)$ is $\varepsilon = \prod_i \alpha_i = \pm 1$.  It is the sign of the functional equation
$$T^{22} P(1/T) = \varepsilon P(T),$$
that is, $a_{22-i} = \varepsilon a_i$ for all $i$.

Since the roots of $P$ that are not $\pm 1$ come in complex conjugate pairs with product $1$, we have $\varepsilon = (-1)^{m_{-1}}$, where $m_{\pm 1}$ is the multiplicity of $\pm 1$ as a root.  As $m_1 + m_{-1}$ is even, this equals $(-1)^{m_1} = (-1)^{\rho}$, where $\rho$ is the {{ KNOWL('k3.fq.picard_rank', 'Picard rank') }}.  So $\varepsilon = -1$ forces $\rho$ to be odd, and in that case $(1-T)(1+T)$ divides $P(T)$.

## k3.fq.real_polynomial

**Title:** Real polynomial of a motive of K3 type

The **real polynomial** of a {{ KNOWL('k3.fq.k3_polynomial', 'motive of K3 type') }} is the polynomial $h(x) \in \Z[x]$ with leading coefficient $q$ whose roots are the real numbers $\alpha_j + \bar\alpha_j = 2\cos(\arg \alpha_j)$, where the $\alpha_j$ run over the normalized Frobenius eigenvalues other than those removed below.  It is defined by
$$P(T) = T^{11} h(T + T^{-1})$$
when the {{ KNOWL('k3.fq.sign', 'sign') }} is $+1$ (so $h$ has degree $11$), and by
$$P(T) = (1 - T^2)\, T^{10} h(T + T^{-1})$$
when the sign is $-1$ (so $h$ has degree $10$).  All roots of $h$ lie in the interval $[-2, 2]$.

## k3.fq.frobenius_angles

**Title:** Frobenius angles of a motive of K3 type

Write the normalized Frobenius eigenvalues of a {{ KNOWL('k3.fq.k3_polynomial', 'motive of K3 type') }} as $\alpha_j = e^{\pi i \theta_j}$.  The **Frobenius angles** are the numbers $\theta_j \in [0, 1]$ for the eigenvalues in the closed upper half plane, listed with multiplicity.  Each pair of complex conjugate eigenvalues $e^{\pm \pi i\theta}$ contributes one angle (displayed as $\pm\theta$), while a root $1$ or $-1$ of multiplicity $m$ contributes $m$ copies of $0$ or $1$ respectively.

The angles are computed numerically; they are displayed to 12 decimal places.

## k3.fq.angle_rank

**Title:** Angle rank of a motive of K3 type

The **angle rank** of a {{ KNOWL('k3.fq.k3_polynomial', 'motive of K3 type') }} with normalized Frobenius eigenvalues $\alpha_i = e^{2\pi i\theta_i}$ is the dimension of the $\Q$-vector space spanned by $\{\theta_1, \ldots, \theta_{22}, 1\}$, minus $1$.  Equivalently, it is the rank of the multiplicative group generated by the $\alpha_i$ (the rank of the free part of this group), and the dimension of the Frobenius torus.

The eigenvalues that are roots of unity contribute nothing, and the remaining $22 - \bar\rho$ eigenvalues (where $\bar\rho$ is the {{ KNOWL('k3.fq.geometric_picard_rank', 'geometric Picard rank') }}) come in pairs $\alpha, \alpha^{-1} = \bar\alpha$.  So the angle rank lies between $0$ and $(22 - \bar\rho)/2$, and it is $0$ exactly when the class is {{ KNOWL('k3.fq.newton_polygon', 'supersingular') }}.

The **angle corank** is $(22 - \bar\rho)/2$ minus the angle rank: the number of independent multiplicative relations among the transcendental eigenvalues, beyond those forced by $\alpha\bar\alpha = 1$.

## k3.fq.newton_polygon

**Title:** Newton polygon of a motive of K3 type

The **Newton polygon** of a {{ KNOWL('k3.fq.k3_polynomial', 'motive of K3 type') }} over $\F_q$ is the $q$-adic Newton polygon of its {{ KNOWL('k3.fq.l_polynomial', 'L-polynomial') }} $L(T)$, that is, the lower convex hull of the points $(i, v_q(c_i))$, where $c_i$ is the coefficient of $T^i$ in $L(T)$ and $v_q = v_p/r$ for $q = p^r$.  Its **slopes** are the $q$-adic valuations of the eigenvalues of Frobenius on $H^2$, and lie in $[0, 2]$.

For a K3 surface the slopes are
$$1 - \tfrac{1}{h}\ (h \text{ times}),\quad 1\ (22 - 2h \text{ times}),\quad 1 + \tfrac{1}{h}\ (h \text{ times}),$$
where $h$ is the {{ KNOWL('k3.fq.height', 'height') }}, or all equal to $1$.  The class has a **K3 Newton polygon** if its slopes have this shape for some $1 \le h \le 11$ or are all equal to $1$.  This is automatic when $q$ is prime, but over non-prime fields there are polynomials whose Newton polygon does not have this shape (for example $P = 2P'$ over $\F_4$, with $P'$ of K3 type over $\F_2$).

A class is
- **ordinary** if $h = 1$, equivalently if $p$ does not divide $a_1$;
- **supersingular** if all slopes equal $1$ (height $h = \infty$), equivalently if every root of $P(T)$ is a root of unity, equivalently if the {{ KNOWL('k3.fq.geometric_picard_rank', 'geometric Picard rank') }} is $22$.

## k3.fq.height

**Title:** Height of a motive of K3 type

The **height** $h$ of a {{ KNOWL('k3.fq.k3_polynomial', 'motive of K3 type') }} is the number of slopes of its {{ KNOWL('k3.fq.newton_polygon', 'Newton polygon') }} that are less than $1$, counted with multiplicity, and $h = \infty$ if all slopes are equal to $1$ (the supersingular case).

For a K3 surface $X$ this is the height of the formal Brauer group $\widehat{\mathrm{Br}}_X$ (Artin–Mazur), and $1 \le h \le 10$ or $h = \infty$.  The surface is ordinary when $h = 1$.  When $h$ is finite, the {{ KNOWL('k3.fq.geometric_picard_rank', 'geometric Picard rank') }} is at most $22 - 2h$.

## k3.fq.newton_elevation

**Title:** Newton elevation of a motive of K3 type

The **Newton elevation** of a {{ KNOWL('k3.fq.k3_polynomial', 'motive of K3 type') }} measures how far its {{ KNOWL('k3.fq.newton_polygon', 'Newton polygon') }} is from the ordinary one, as for {{ KNOWL('av.fq.newton_elevation', 'abelian varieties') }}.  It is the number of lattice points $(x, y)$ with $0 \le x \le 11$ that lie on or above the ordinary Newton polygon (with vertices $(0,0)$, $(1,0)$, $(21,20)$, $(22,22)$) and strictly below the Newton polygon of $L(T)$.

The **Newton coelevation** is the number of lattice points $(x, y)$ with $0 \le x \le 11$ that lie on or above the Newton polygon of $L(T)$ and strictly below the supersingular polygon $y = x$.

For a {{ KNOWL('k3.fq.newton_polygon', 'K3 Newton polygon') }} of {{ KNOWL('k3.fq.height', 'height') }} $h$, the elevation is $h - 1$ and the coelevation is $12 - h$; in the supersingular case they are $11$ and $0$.  Both are left undefined when the Newton polygon is not a K3 Newton polygon (possible only over non-prime fields).

## k3.fq.picard_rank

**Title:** Picard rank of a motive of K3 type

The **Picard rank** $\rho$ of a {{ KNOWL('k3.fq.k3_polynomial', 'motive of K3 type') }} with polynomial $P(T)$ is the multiplicity of $T = 1$ as a root of $P(T)$.

For a K3 surface $X/\F_q$, the Tate conjecture (a theorem for K3 surfaces over finite fields) says that $\rho$ is the rank of the Néron–Severi group $\mathrm{NS}(X)$.  Over the extension $\F_{q^n}$ the Picard rank is $\rho_n$, the multiplicity of $1$ as a root of the base change $P_n(T) = q^n\prod(1 - \alpha_i^n T)$; it only changes at the {{ KNOWL('k3.fq.geometric_extension_degree', 'extension degrees') }} determined by the {{ KNOWL('k3.fq.cyclotomic_part', 'cyclotomic part') }} of $P$, and is eventually equal to the {{ KNOWL('k3.fq.geometric_picard_rank', 'geometric Picard rank') }}.

The parity of $\rho$ is determined by the {{ KNOWL('k3.fq.sign', 'sign') }}: $\varepsilon = (-1)^{\rho}$.  A projective surface has $\rho \ge 1$.

## k3.fq.geometric_picard_rank

**Title:** Geometric Picard rank of a motive of K3 type

The **geometric Picard rank** $\bar\rho$ of a {{ KNOWL('k3.fq.k3_polynomial', 'motive of K3 type') }} is the number of roots of $P(T)$ (with multiplicity) that are roots of unity, that is, the degree of its {{ KNOWL('k3.fq.cyclotomic_part', 'cyclotomic part') }} $C(T)$.

For a K3 surface $X/\F_q$ it is the rank of $\mathrm{NS}(X_{\overline{\F}_q})$, by the Tate conjecture.  It is always even, since the other roots come in complex conjugate pairs; it equals $22$ exactly when the class is {{ KNOWL('k3.fq.newton_polygon', 'supersingular') }}, and otherwise it is at most $22 - 2h$, where $h$ is the {{ KNOWL('k3.fq.height', 'height') }}.

## k3.fq.cyclotomic_part

**Title:** Cyclotomic part of a motive of K3 type

The **cyclotomic part** $C(T)$ of the polynomial $P(T)$ of a {{ KNOWL('k3.fq.k3_polynomial', 'motive of K3 type') }} is the product of its cyclotomic factors, normalized so that $C(0) = 1$:
$$C(T) = \pm\prod_m \Phi_m(T)^{e_m},$$
where $\Phi_m$ is the $m$th cyclotomic polynomial (so $\Phi_1(T) = T - 1$ and $\Phi_2(T) = T + 1$).  The LMFDB displays the factorization as $\prod_m \Phi_m^{e_m}$.  By the Tate conjecture, the roots of $C$ correspond to the algebraic classes in $H^2$; the multiplicity $e_1$ of $\Phi_1$ is the {{ KNOWL('k3.fq.picard_rank', 'Picard rank') }} and the degree of $C$ is the {{ KNOWL('k3.fq.geometric_picard_rank', 'geometric Picard rank') }}.

The quotient $P_T = P/C$ is the {{ KNOWL('k3.fq.transcendental_part', 'transcendental part') }}.

## k3.fq.geometric_extension_degree

**Title:** Geometric extension degree of a motive of K3 type

The **geometric extension degree** $E$ of a {{ KNOWL('k3.fq.k3_polynomial', 'motive of K3 type') }} is the least common multiple of the $m$ such that the cyclotomic polynomial $\Phi_m$ divides $P(T)$ (and $E = 1$ if there are none).  It is the smallest $n$ such that the {{ KNOWL('k3.fq.picard_rank', 'Picard rank') }} $\rho_n$ over $\F_{q^n}$ equals the {{ KNOWL('k3.fq.geometric_picard_rank', 'geometric Picard rank') }}; for a K3 surface, all of $\mathrm{NS}(X_{\overline{\F}_q})$ is defined over $\F_{q^E}$.

More precisely, if $e_m$ is the multiplicity of $\Phi_m$, then
$$\rho_n = \sum_{m \mid n} e_m\, \varphi(m).$$
So $\rho_n$ only depends on $\mathrm{lcm}\{m : m \mid n,\ e_m > 0\}$, and the **extension degrees** at which the Picard rank can change are the least common multiples of subsets of $\{m : e_m > 0\}$ (including $1$).

## k3.fq.artin_tate

**Title:** Artin–Tate formula for motives of K3 type

Let $X$ be a K3 surface over $\F_q$ with {{ KNOWL('k3.fq.picard_rank', 'Picard rank') }} $\rho$ and polynomial $P(T)$.  The Artin–Tate formula (a theorem for K3 surfaces, since the Tate conjecture holds for them) states that
$$\left.\frac{P(T)}{(1-T)^{\rho}}\right|_{T=1} = \#\mathrm{Br}(X)\cdot\lvert\operatorname{disc}\mathrm{NS}(X)\rvert,$$
where $\mathrm{Br}(X)$ is the Brauer group.

For any {{ KNOWL('k3.fq.k3_polynomial', 'motive of K3 type') }} and $n \ge 1$, the LMFDB records the **Artin–Tate value**
$$AT_n = \left.\frac{P_n(T)}{(1-T)^{\rho_n}}\right|_{T=1},$$
where $P_n(T) = q^n \prod_i (1 - \alpha_i^n T)$ is the polynomial of the base change to $\F_{q^n}$ and $\rho_n$ its Picard rank, for each of the {{ KNOWL('k3.fq.geometric_extension_degree', 'extension degrees') }} $n$ at which the Picard rank can change.  It is a positive integer.

Since $\#\mathrm{Br}(X)$ is a perfect square, the squarefree part of $AT_1$ gives the {{ KNOWL('k3.fq.ns_disc_class', 'discriminant class') }} of $\mathrm{NS}(X)$.

## k3.fq.ns_disc_class

**Title:** Néron–Severi discriminant class

The **NS discriminant class** of a {{ KNOWL('k3.fq.k3_polynomial', 'motive of K3 type') }} with {{ KNOWL('k3.fq.picard_rank', 'Picard rank') }} $\rho \ge 1$ is the squarefree integer
$$(-1)^{\rho - 1}\, \mathrm{sqf}(AT_1),$$
where $AT_1$ is the {{ KNOWL('k3.fq.artin_tate', 'Artin–Tate value') }} and $\mathrm{sqf}$ denotes the squarefree part.

For a K3 surface $X$, the Artin–Tate formula and the fact that $\#\mathrm{Br}(X)$ is a square show that this is the class of $\operatorname{disc}\mathrm{NS}(X)$ in $\Q^\times/\Q^{\times 2}$ (the sign comes from the signature $(1, \rho - 1)$ of $\mathrm{NS}(X)$).  It is not defined when $\rho = 0$.

## k3.fq.ej_condition

**Title:** Artin–Tate condition of Elsenhans and Jahnel

Let $P(T)$ be the polynomial of a {{ KNOWL('k3.fq.k3_polynomial', 'motive of K3 type') }} with {{ KNOWL('k3.fq.picard_rank', 'Picard rank') }} $\rho$, and let $S(T) = P(T)/(1-T)^{\rho}$.  The class satisfies the **Artin–Tate condition** of Elsenhans and Jahnel if $S(-1)$ is a perfect square (where $0$ counts as a square), as tested by Kedlaya and Sutherland \cite{mr:3540942,arxiv:1511.06945} (`k3-scripts/ej-test.sage` in root-unitary).

When $1 + T$ does not divide $P(T)$, the condition says that the product $AT_1 AT_2$ of the {{ KNOWL('k3.fq.artin_tate', 'Artin–Tate values') }} over $\F_q$ and $\F_{q^2}$ is a square.  Elsenhans and Jahnel showed that the polynomial of a K3 surface satisfies this condition, so it is a necessary condition for the class to contain a K3 surface.  [Add the Elsenhans–Jahnel reference.]

## k3.fq.transcendental_part

**Title:** Transcendental part of a motive of K3 type

The **transcendental part** of the polynomial $P(T)$ of a {{ KNOWL('k3.fq.k3_polynomial', 'motive of K3 type') }} is
$$P_T(T) = P(T)/C(T) \in \Z[T],$$
where $C(T)$ is the {{ KNOWL('k3.fq.cyclotomic_part', 'cyclotomic part') }}.  It has constant term $q$, degree $22 - \bar\rho$ (where $\bar\rho$ is the {{ KNOWL('k3.fq.geometric_picard_rank', 'geometric Picard rank') }}), and none of its roots is a root of unity.  For a K3 surface it is the polynomial of Frobenius acting on the transcendental part of $H^2$, the orthogonal complement of the algebraic classes.  For a {{ KNOWL('k3.fq.newton_polygon', 'supersingular') }} class it is the constant $q$.

The LMFDB lists the factorization of $P_T$ over $\Q$ into primitive irreducible factors with positive constant term, with multiplicities.  When $q$ is prime, $P_T$ is always irreducible: since the constant terms of the factors multiply to (a divisor of) $q$, all but one factor would have constant term $1$, and by Kronecker's theorem such a factor would be a product of cyclotomic polynomials.  Over non-prime fields $P_T$ can be a power of an irreducible polynomial (for instance for a base change) or a product of distinct irreducible polynomials.

## k3.fq.transcendental_number_field

**Title:** Number field of the transcendental part

For each irreducible factor $f$ of the {{ KNOWL('k3.fq.transcendental_part', 'transcendental part') }} of a {{ KNOWL('k3.fq.k3_polynomial', 'motive of K3 type') }}, the LMFDB records the number field $K = \Q[T]/(f)$.  Since the roots of $f$ lie on the unit circle and are not real, $K$ is a {{ KNOWL('nf.cm_field', 'CM field') }}, with maximal totally real subfield generated by $\alpha + \alpha^{-1}$.

The field is linked to its page in the LMFDB when it is in the {{ KNOWL('nf', 'number field database') }}; otherwise it is described by a reduced defining polynomial (computed with `polredabs`, or `polredbest` when `polredabs` does not finish within the time limit).

## k3.fq.transcendental_galois_group

**Title:** Galois group of the transcendental part

For each irreducible factor $f$ of the {{ KNOWL('k3.fq.transcendental_part', 'transcendental part') }} of a {{ KNOWL('k3.fq.k3_polynomial', 'motive of K3 type') }}, this is the Galois group of $f$ over $\Q$, given as a {{ KNOWL('gg.label', 'transitive permutation group') }} $n$T$k$ where $n = \deg f$.

When searching by Galois group, a class matches if the Galois group of one of the factors of its transcendental part is among the groups given.  A group can be given by its transitive group label (such as 4T3) or by a {{ KNOWL('nf.galois_group.name', 'name') }} (such as D4).

## k3.fq.point_counts

**Title:** Point counts of varieties of K3 type

Let $\alpha_1, \ldots, \alpha_{22}$ be the normalized Frobenius eigenvalues of a {{ KNOWL('k3.fq.k3_polynomial', 'motive of K3 type') }} $M$ over $\F_q$.  By the Lefschetz trace formula, a K3 surface $X$ with $H^2(X)(1)$ in this class has
$$\#X(\F_{q^n}) = 1 + q^{2n} + q^n \sum_{i} \alpha_i^n.$$
In particular $\#X(\F_q) = 1 + q^2 - a_1$.

The LMFDB also lists the point counts of fourfolds of K3 type whose K3-type piece is in the class, assuming that the rest of their cohomology consists of Tate classes defined over $\F_q$:

- a cubic fourfold $Y$ ($b_2 = 1$, $b_4 = 23$, $M = H^4_{\mathrm{prim}}(Y)(1)$):
$$1 + q^n + q^{2n} + q^{3n} + q^{4n} + q^{2n}\sum_i\alpha_i^n;$$
- a Gushel–Mukai fourfold ($b_2 = 1$, $b_4 = 24$, $M = H^4_{\mathrm{van}}(Y)(1)$):
$$1 + q^n + 2q^{2n} + q^{3n} + q^{4n} + q^{2n}\sum_i\alpha_i^n;$$
- a Verra fourfold, a double cover of $\mathbb{P}^2 \times \mathbb{P}^2$ branched along a divisor of bidegree $(2,2)$ ($\chi = 30$, $b_2 = b_6 = 2$, $b_4 = 24$):
$$1 + 2q^n + 2q^{2n} + 2q^{3n} + q^{4n} + q^{2n}\sum_i\alpha_i^n.$$
Here the 22-dimensional K3-type piece of $H^4(Y)(1)$ is taken to contain one of the three algebraic classes pulled back from $\mathbb{P}^2\times\mathbb{P}^2$.  [Check against the intended definition of the K3-type piece of a Verra fourfold.]

These are virtual counts: for classes that do not contain such a variety they may be negative.

## k3.fq.has_nonnegative_counts

**Title:** Nonnegative point counts

A {{ KNOWL('k3.fq.k3_polynomial', 'motive of K3 type') }} over $\F_q$ has **nonnegative counts** if the {{ KNOWL('k3.fq.point_counts', 'point counts') }} $1 + q^{2n} + q^n\sum_i\alpha_i^n$ of a K3 surface in the class are nonnegative for all $n \ge 1$.  Since $|\sum_i \alpha_i^n| \le 22$, only the finitely many $n$ with $q^n < 22$ need to be checked.

This is a necessary condition for the class to contain a K3 surface.

## k3.fq.has_count_growth

**Title:** Point count growth

A {{ KNOWL('k3.fq.k3_polynomial', 'motive of K3 type') }} over $\F_q$ has **count growth** if, for a K3 surface $X$ in the class, the number of closed points of degree $n$,
$$\frac{1}{n}\sum_{d \mid n}\mu(n/d)\,\#X(\F_{q^d}),$$
is nonnegative for every $n \ge 1$.  This is a necessary condition for the class to contain a K3 surface.  It implies that $\#X(\F_{q^m}) \ge \#X(\F_{q^n})$ whenever $n \mid m$ (for example, a K3 surface has at least as many points over $\F_{16}$ as over $\F_4$), and in particular that the {{ KNOWL('k3.fq.has_nonnegative_counts', 'point counts are nonnegative') }}.

## k3.fq.pointless

**Title:** Pointless motive of K3 type

A {{ KNOWL('k3.fq.k3_polynomial', 'motive of K3 type') }} over $\F_q$ is **pointless** if a K3 surface in the class has no $\F_q$-points, that is, if
$$\#X(\F_q) = 1 + q^2 - a_1 = 0.$$
[Add references: pointless K3 surfaces over $\F_2$; the enumeration of their Néron–Severi lattices by Brandhorst.]

## k3.fq.realizability

**Title:** Realizability of a motive of K3 type

A class of {{ KNOWL('k3.fq.k3_polynomial', 'motives of K3 type') }} over $\F_q$ is **realizable** by a given type of object if there is such an object over $\F_q$ whose K3-type motive has the polynomial $P(T)$ of the class.  The LMFDB records, for each type, whether this is known to hold (*yes*), known to fail (*no*), or *unknown*; it is *not computed* if the question has not been considered yet.

The types are:
- **K3 surfaces** $X$, with $M = H^2(X)(1)$;
- **cubic fourfolds** $Y$, with $M = H^4_{\mathrm{prim}}(Y)(1)$.  A realization comes with a positive definite marking: an embedding of a lattice with Gram matrix $\begin{pmatrix} 3 & a \\ a & b\end{pmatrix}$, of discriminant $3b - a^2$, into the algebraic classes (compare Hassett's special cubic fourfolds);
- **Gushel–Mukai fourfolds**, with $M = H^4_{\mathrm{van}}(Y)(1)$;
- **Verra fourfolds**, with $M$ the K3-type piece of $H^4(Y)(1)$ (see {{ KNOWL('k3.fq.point_counts', 'point counts') }});
- **Kummer surfaces** $\mathrm{Km}(A)$ of abelian surfaces $A/\F_q$ (see {{ KNOWL('k3.fq.kummer', 'Kummer candidates') }});
- the **noncommutative** versions of the first four: K3 categories over $\F_q$, such as twisted K3 surfaces and the Kuznetsov components of cubic, Gushel–Mukai and Verra fourfolds, with $M$ the K3-type part of the Mukai motive.

Necessary conditions for a class to contain a K3 surface include $\rho \ge 1$ (where $\rho$ is the {{ KNOWL('k3.fq.picard_rank', 'Picard rank') }}), the {{ KNOWL('k3.fq.ej_condition', 'Artin–Tate condition') }} and {{ KNOWL('k3.fq.has_count_growth', 'count growth') }}.  For ordinary classes the question can be decided using Taelman's description of ordinary K3 surfaces over finite fields in terms of their canonical lifts.  [Describe the methods used for the other types, and add references.]

## k3.fq.kummer

**Title:** Kummer candidates

For $p > 2$, the Kummer surface $\mathrm{Km}(A)$ of an abelian surface $A/\F_q$ has
$$P(T) = P_{\wedge^2}(T) \cdot P_{16}(T),$$
where $P_{\wedge^2}$ is the polynomial of the normalized $\wedge^2 H^1(A)(1)$ (whose roots are the $\pi_i\pi_j/q$ for the Frobenius eigenvalues $\pi_1, \ldots, \pi_4$ of $A$), and $P_{16}$ is the characteristic polynomial of the permutation action of Frobenius on the 16 exceptional curves, which are indexed by the points of $A[2]$.

The **Kummer candidates** of a {{ KNOWL('k3.fq.k3_polynomial', 'motive of K3 type') }} are the {{ KNOWL('av.isogeny_class', 'isogeny classes') }} of abelian surfaces $A/\F_q$ whose normalized wedge square has the same {{ KNOWL('k3.fq.transcendental_part', 'transcendental part') }} as the class: these are the candidates for a Kummer structure.  The list is empty if there are none; it is not computed in characteristic 2 or when the LMFDB has no abelian surfaces over $\F_q$.  A class with no Kummer candidates, or one that fails the necessary conditions for being the class of a K3 surface, is not the class of a Kummer surface.

## k3.fq.polarization_degrees

**Title:** Polarization degrees

For an ordinary {{ KNOWL('k3.fq.k3_polynomial', 'motive of K3 type') }}, the **polarization degrees** are the integers $2d \le 100$ such that the class contains a K3 surface with a primitive polarization of degree $2d$.  They are determined using Taelman's description of ordinary K3 surfaces over finite fields.  They are not computed for other classes.

## k3.fq.primitive

**Title:** Primitive motive of K3 type

A class of {{ KNOWL('k3.fq.k3_polynomial', 'motives of K3 type') }} over $\F_q$ is **primitive** if it is not the base change of a class over a proper subfield of $\F_q$, that is, if there is no class over $\F_{q_0}$, with $q_0^d = q$ and $d > 1$, whose base change
$$P_d(T) = q_0^d\prod_i (1 - \alpha_i^d T)$$
is the polynomial of the class.  For a class that is not primitive, the LMFDB lists the classes over subfields (its **primitive models**) whose base change it is.

## k3.fq.twist

**Title:** Twists of motives of K3 type

Two classes of {{ KNOWL('k3.fq.k3_polynomial', 'motives of K3 type') }} over $\F_q$ are **twists** of each other if they become equal over some finite extension $\F_{q^n}$, that is, if their base changes $P_n$ and $P'_n$ agree.  This is an equivalence relation, and its classes are called **twist classes**; all members of a twist class have the same $q$.

The representative of a twist class (the **geometric class**) is its member with the smallest label.  For each class the LMFDB lists the number of classes in its twist class (including the class itself) and, when this number is at most $50$, the list of twists: for each twist, a minimal degree $n$ (for divisibility) over which the two classes become equal, and the label of their common base change to $\F_{q^n}$, as for {{ KNOWL('av.twist', 'twists of abelian varieties') }}.  Supersingular classes and classes with large {{ KNOWL('k3.fq.cyclotomic_part', 'cyclotomic part') }} have very large twist classes; use a search on the geometric class to list them.  The **maximum twist degree** is the largest $n$ in the list.

The geometric class is computed from the multiset of $\{M\theta_j\} \bmod 1$ over the transcendental Frobenius angles $\theta_j = \arg(\alpha_j)/2\pi$, with $M = \mathrm{lcm}\{m : \varphi(m) \le 44\}$; the minimal degrees are then verified by exact base change.

## k3.fq.lattices

**Title:** Candidate Néron–Severi lattices

Let $X/\F_q$ be a K3 surface whose class has polynomial $P(T)$ with {{ KNOWL('k3.fq.cyclotomic_part', 'cyclotomic part') }} $C(T)$ of degree $\bar\rho$.  Then $N = \mathrm{NS}(X_{\overline{\F}_q})$ is an even lattice of signature $(1, \bar\rho - 1)$, and Frobenius acts on it by an isometry $\varphi$ of finite order with characteristic polynomial determined by $C(T)$, preserving the ample cone (which is a Weyl chamber of the positive cone).  The invariant sublattice $N^\varphi$ is $\mathrm{NS}(X)$.

In some cases $P$ determines the pair $(N, \varphi)$ up to isomorphism; in general there are finitely many possibilities.  The LMFDB lists the isomorphism classes of pairs $(N, \varphi)$ that are not ruled out by the following necessary conditions:

- *Artin–Tate*: the orders $AT_n/\lvert\operatorname{disc} N^{\varphi^n}\rvert$ of the Brauer groups (see {{ KNOWL('k3.fq.artin_tate', 'Artin–Tate formula') }}) are perfect squares;
- *$\ell$-adic*: for each prime $\ell \ne p$, $N \otimes \Z_\ell$ embeds primitively in the even unimodular $\Z_\ell$-lattice of rank 22, with complement carrying the transcendental part;
- *crystalline*: compatibility with the Newton polygon (Ogus' conditions on $N \otimes \Z_p$);
- *chamber*: $\varphi$ preserves a Weyl chamber.

Each candidate need not come from an actual K3 surface.  For each candidate the LMFDB records the genus of $N$, the discriminant and rank of $N^\varphi$, the order of $\mathrm{Br}(X)$, the minimal degree $v^2$ of a big and nef class over $\F_q$ (the minimal positive norm in $N^\varphi$) and of an ample class (a class in $N^\varphi$ in the interior of a $\varphi$-stable Weyl chamber), whether $N$ contains roots (vectors of norm $-2$, corresponding to smooth rational curves), and whether $N$ primitively contains the Kummer lattice (by Nikulin, the lattice-theoretic condition for $X_{\overline{\F}_q}$ to be a Kummer surface).

The candidates were enumerated with Oscar's machinery for lattices with isometry (Brandhorst–Muller).  The computation is expensive, and it has only been carried out for some classes.

## k3.fq.search_input

**Title:** Find a motive of K3 type

To find a specific class of {{ KNOWL('k3.fq.k3_polynomial', 'motives of K3 type') }} over a finite field, enter either its {{ KNOWL('k3.fq.lmfdb_label', 'label') }}, such as `2.p.abs_ru_aeom_vqs_aczxs_imtm_atkou_bkkbs_aceprc_cvpea_adchca`, or its polynomial: either $P(T)$ (with constant term $q$) or the {{ KNOWL('k3.fq.l_polynomial', 'L-polynomial') }} $L(T) = q^{-1}P(qT)$ (with constant term $1$), as a polynomial of degree $22$ in one variable, such as `2*(1-T)^22`, or as the list of its 23 coefficients, constant term first.

## k3.fq.initial_coefficients

**Title:** Initial coefficients of a motive of K3 type

The LMFDB allows one to search for {{ KNOWL('k3.fq.k3_polynomial', 'motives of K3 type') }} over finite fields by specifying the first few coefficients $a_1, a_2, \ldots$ of the polynomial
$$P(T) = q + a_1 T + a_2 T^2 + \cdots + a_{22}T^{22},$$
in increasing order of degree, as a list such as `[-1, 3]`.  The constant term $a_0 = q$ is **not** included.  Note that these are the coefficients of $P(T)$, as in the {{ KNOWL('k3.fq.lmfdb_label', 'label') }}, rather than those of the {{ KNOWL('k3.fq.l_polynomial', 'L-polynomial') }}, whose coefficient of $T^i$ is $q^{i-1}a_i$.

The search is much faster when the base field is specified as well.

## rcs.source.k3.fq

**Title:** Source of the data on motives of K3 type over finite fields

The polynomials were enumerated using Kedlaya's algorithm for enumerating Weil polynomials (as implemented in root-unitary and in SageMath's `WeilPolynomials`), as in the work of Kedlaya and Sutherland on zeta functions of K3 surfaces \cite{mr:3540942,arxiv:1511.06945}.  Polynomials with sign $-1$ were enumerated as $(1 - T^2)P'(T)$, with $P'$ of degree $20$ and sign $+1$.

The invariants were computed in SageMath and PARI/GP.  [Describe the computation of twists, Kummer candidates and realizability.]  The candidate Néron–Severi lattices were enumerated using Oscar's machinery for lattices with isometry, due to Brandhorst and Muller, [with Magma used for ...].

[Add the authors of the computation, the date, and links to the code.]

## rcs.ack.k3.fq

**Title:** Acknowledgments for motives of K3 type over finite fields

[To be written: who computed the data, and funding.]

The computations build on the work of Kedlaya and Sutherland \cite{mr:3540942,arxiv:1511.06945} and on Oscar's machinery for lattices with isometry by Brandhorst and Muller.

## rcs.cite.k3.fq

**Title:** How to cite motives of K3 type over finite fields

To cite the data on motives of K3 type over finite fields, please cite the LMFDB, and [the paper describing the computation, once available].

## rcs.cande.k3.fq

**Title:** Completeness of the data on motives of K3 type over finite fields

For each base field $\F_q$ represented in the database, the database contains every polynomial $P(T) = q + a_1 T + \cdots \pm q T^{22} \in \Z[T]$ with $a_{22-i} = \pm a_i$ and all roots on the unit circle (both signs), regardless of whether it is the polynomial of an actual K3 surface.  The base fields currently included are [$q = 2$, ...].

Some quantities have not been computed for every class, and are shown as *not computed*:
- the LMFDB labels of the {{ KNOWL('k3.fq.transcendental_number_field', 'number fields') }} of the transcendental part (pending a lookup in the number field database);
- {{ KNOWL('k3.fq.realizability', 'realizability') }}, beyond the necessary conditions;
- the {{ KNOWL('k3.fq.polarization_degrees', 'polarization degrees') }};
- the list of {{ KNOWL('k3.fq.twist', 'twists') }} for classes with more than $50$ twists;
- the {{ KNOWL('k3.fq.lattices', 'candidate Néron–Severi lattices') }}, which have only been enumerated for [...].

## rcs.rigor.k3.fq

**Title:** Reliability of the data on motives of K3 type over finite fields

The enumeration of the polynomials and the computation of their factorizations, cyclotomic and transcendental parts, Picard ranks, Artin–Tate values, Newton polygons, point counts and base changes use exact arithmetic and do not depend on any unproven assumptions.

- The {{ KNOWL('k3.fq.frobenius_angles', 'Frobenius angles') }} are computed numerically.  [State how the {{ KNOWL('k3.fq.angle_rank', 'angle rank') }} is computed, e.g. numerically with PARI's `lindep` at a given precision.]
- The {{ KNOWL('k3.fq.twist', 'geometric class') }} is computed from numerical Frobenius angles, and the twists are then verified by exact base change.
- The identification of the Picard rank with the rank of the Néron–Severi group, and the Artin–Tate formula, rely on the Tate conjecture, which is a theorem for K3 surfaces over finite fields.
- The {{ KNOWL('k3.fq.lattices', 'candidate lattices') }} are only required to satisfy the listed necessary conditions; they need not be realized by K3 surfaces.  [Describe the cross-checks between Magma and Oscar, and the known bugs that were patched.]
