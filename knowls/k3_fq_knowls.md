# Remaining draft knowls for motives of K3 type over finite fields

The knowls used by the search form and the home page of the `k3fq` section of
the LMFDB (`lmfdb/k3fq/` on the `k3fq` branch, url `/Motive/K3/Fq/`) are in
`k3_knowls.yaml` and `k3_rcs_knowls.yaml`, in the format read by
`db.knowldb.yaml_import`.  This file has the remaining drafts, used only on the
pages of individual classes; they have not been converted or shortened yet.  Each section gives
the knowl id, a proposed title, and the proposed content (LMFDB knowl markdown,
with `{{ KNOWL(...) }}` cross references and `\cite{...}` as in existing knowls).
Text in [square brackets] marks places that need input.  Knowls for interesting
classes follow the usual pattern `k3.fq.<label>.top`, and are listed on
`/Motive/K3/Fq/interesting`.

Suggested edit to an existing knowl: the page `/Motive` renders the knowl
`motives.about`, which currently lists only hypergeometric motives.  Add

    <li><a href="{{ url_for('k3fq.index') }}">Motives of K3 type over finite fields</a>

after the hypergeometric line.

---

## k3.fq.real_polynomial

**Title:** Real polynomial of a motive of K3 type

The **real polynomial** of a {{ KNOWL('k3.fq.motive', 'motive of K3 type') }} is the polynomial $h(x) \in \Z[x]$ with leading coefficient $q$ whose roots are the real numbers $\alpha_j + \bar\alpha_j = 2\cos(\arg \alpha_j)$, where the $\alpha_j$ run over the normalized Frobenius eigenvalues other than those removed below.  It is defined by
$$P(T) = T^{11} h(T + T^{-1})$$
when the {{ KNOWL('k3.fq.sign', 'sign') }} is $+1$ (so $h$ has degree $11$), and by
$$P(T) = (1 - T^2)\, T^{10} h(T + T^{-1})$$
when the sign is $-1$ (so $h$ has degree $10$).  All roots of $h$ lie in the interval $[-2, 2]$.

## k3.fq.frobenius_angles

**Title:** Frobenius angles of a motive of K3 type

Write the normalized Frobenius eigenvalues of a {{ KNOWL('k3.fq.motive', 'motive of K3 type') }} as $\alpha_j = e^{\pi i \theta_j}$.  The **Frobenius angles** are the numbers $\theta_j \in [0, 1]$ for the eigenvalues in the closed upper half plane, listed with multiplicity.  Each pair of complex conjugate eigenvalues $e^{\pm \pi i\theta}$ contributes one angle (displayed as $\pm\theta$), while a root $1$ or $-1$ of multiplicity $m$ contributes $m$ copies of $0$ or $1$ respectively.

The angles are computed numerically; they are displayed to 12 decimal places.

## k3.fq.kummer

**Title:** Kummer candidates

For $p > 2$, the Kummer surface $\mathrm{Km}(A)$ of an abelian surface $A/\F_q$ has
$$P(T) = P_{\wedge^2}(T) \cdot P_{16}(T),$$
where $P_{\wedge^2}$ is the polynomial of the normalized $\wedge^2 H^1(A)(1)$ (whose roots are the $\pi_i\pi_j/q$ for the Frobenius eigenvalues $\pi_1, \ldots, \pi_4$ of $A$), and $P_{16}$ is the characteristic polynomial of the permutation action of Frobenius on the 16 exceptional curves, which are indexed by the points of $A[2]$.

The **Kummer candidates** of a {{ KNOWL('k3.fq.motive', 'motive of K3 type') }} are the {{ KNOWL('av.isogeny_class', 'isogeny classes') }} of abelian surfaces $A/\F_q$ whose normalized wedge square has the same {{ KNOWL('k3.fq.transcendental_part', 'transcendental part') }} as the class: these are the candidates for a Kummer structure.  The list is empty if there are none; it is not computed in characteristic 2 or when the LMFDB has no abelian surfaces over $\F_q$.  A class with no Kummer candidates, or one that fails the necessary conditions for being the class of a K3 surface, is not the class of a Kummer surface.

## k3.fq.polarization_degrees

**Title:** Polarization degrees

For an ordinary {{ KNOWL('k3.fq.motive', 'motive of K3 type') }}, the **polarization degrees** are the integers $2d \le 100$ such that the class contains a K3 surface with a primitive polarization of degree $2d$.  They are determined using Taelman's description of ordinary K3 surfaces over finite fields, subject to the same existence hypotheses as for {{ KNOWL('k3.fq.realizability', 'realizability') }}.  They are not computed for other classes.
