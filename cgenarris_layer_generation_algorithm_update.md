# cgenarris Layer (Film) Generation — Algorithm Update

This document records the design of the new film-generation algorithm for cgenarris.
Section 1 summarizes the generation logic of the code as it currently exists.
Section 2 specifies the new algorithm: the high-level workflow, followed by the
mathematical and algorithmic details of each component. The document is focused on
design, algorithm, and mathematical logic; code is referenced only where necessary.

---

## 1. Current generation logic

The current entry point is `mpi_generate_layer_with_vdw_cutoff_matrix`. It generates
2D-periodic molecular films whose in-plane lattices are integer supercells of a given
substrate surface lattice, organized by layer group, in parallel with MPI.

### 1.1 Setup

- Inputs: a pairwise van-der-Waals cutoff matrix (per atom pair), number of structures
  per layer group, Z (number of molecules per cell), target volume mean/std, interface
  area mean/std, a volume multiplier, a tolerance `TOL`, max attempts, a space-group
  distribution type, the substrate 2D lattice vectors, a volume-resampling cadence, and
  a random seed.
- Each MPI rank reads the molecule from `geometry.in`, recenters it to its centroid,
  and seeds its own random number generators.
- An initial target volume is drawn from a Gaussian `N(volume_mean, volume_std)`.

### 1.2 Substrate supercell ("epitaxy") combinations

All integer 2×2 transformation matrices `[[f1, f2], [f3, f4]]` applied to the substrate
lattice vectors are enumerated by brute force over a 4-dimensional integer range. A
combination is kept if it is right-handed (positive determinant), its cell area is below
`0.5 × target_volume`, the aspect ratio |a|/|b| is below 3, and the cell angle is within
[30°, 150°]. The surviving integer 4-tuples are stored in a preallocated array of
100,000,000 × 4 ints (~1.6 GB per rank).

### 1.3 Compatible layer group detection (trial-and-error)

The molecule's internal symmetry is never determined explicitly. Instead:

1. The atoms of the same element lying on the farthest distance shell from the centroid
   are collected ("equivalent atoms").
2. Candidate molecular axes are constructed from pairs of those atoms (bisectors and
   cross products), deduplicated.
3. For each of the 80 layer groups, a **fake lattice** with large volume is constructed.
   For each Wyckoff position with multiplicity equal to Z, the molecule is placed at the
   position, rotated through every (candidate axis × 16 hardcoded viewing directions)
   combination, all layer-group operations are applied, and the position is declared
   compatible if the symmetry-generated copies overlap the original within `TOL` and
   merge to the correct multiplicity. The working (axis, viewing direction) pairs are
   cached and reused at generation time.

A by-product of this numeric procedure is the **overlap list** of each compatible
special position: the indices of the symmetry-generated copies that coincide with the
first molecule, used later to verify placements and merge coincident copies.

### 1.4 Per-layer-group generation loop

For each compatible layer group, until the per-group quota is met or `max_attempts` is
exhausted, every rank repeatedly:

1. Samples substrate combinations at random until one matches the layer group's 2D
   lattice class (oblique / rectangular / square / hexagonal, tested with tolerances on
   |a| = |b| and γ = 90°/120°) and its area lies within `interface_area_mean ± std`.
2. Builds the out-of-plane vector c so the cell volume equals
   `target_volume × volume_multiplier`. The c direction follows the 3D Bravais class of
   the layer group: continuous tilt (two random angles) for the triclinic/oblique class,
   tilt within the bc-plane for the rectangular/monoclinic class, perpendicular c for
   the remaining classes.
3. Picks one compatible Wyckoff position at random, places the molecule with a random
   position along the position's free coordinates and a random orientation drawn from the
   cached (axis, viewing direction) pairs, and applies all layer-group operations.
4. Checks all interatomic distances (including periodic images) against the vdW cutoff
   matrix; on success the structure is accepted.
5. Periodically resamples the target volume from the Gaussian.

### 1.5 MPI collection

After each batch of attempts, all ranks gather their success flags to rank 0; rank 0
receives successful structures, writes them to `geometry.out` (lattice, coordinates,
layer group, Wyckoff position, COM, Euler angles, epitaxy matrix), increments the shared
counter, and broadcasts stop/success flags. All ranks synchronize at layer-group
boundaries.

### 1.6 Limitations motivating the redesign

1. **Implicit, tolerance-fragile symmetry handling.** Compatibility depends on numeric
   overlap tests on a fake lattice; the molecular point group is never known explicitly,
   and results depend on `TOL` and on 16 hardcoded viewing directions.
2. **Redundant epitaxy enumeration.** Infinitely many integer matrices describe the same
   sublattice (unimodular redundancy); the brute-force table stores duplicates and costs
   ~1.6 GB per rank, and random sampling from it is silently biased toward sublattices
   with more representations.
3. **Volume-first lattice construction.** c is stretched to meet a sampled volume
   (`volume_multiplier` hack); films can contain unphysical vacuum along z, which is
   meaningless for future multilayer stacking and wastes attempts on z-clashing or
   z-padded structures.
4. **Tight MPI coupling.** Gather/broadcast every batch; ranks idle at synchronization
   points; early-stop logic complicates the control flow.
5. **Rejection-driven distribution shaping.** The volume distribution is imposed by
   rejection and resampling cadence (`vol_attempt`) rather than by construction.

---

## 2. New algorithm

### 2.0 Design goals

- Determine molecular symmetry **explicitly and exactly** (up to a distance tolerance),
  and use it symbolically: no fake lattices, no trial orientations.
- Enumerate epitaxy matrices **canonically** (each distinct sublattice exactly once),
  with physically motivated area bounds.
- Construct the out-of-plane vector **z-first**: the film must be correctly packed along
  z (contact condition) *before* any volume consideration, because cells will later be
  replicated along c for multilayer generation. Volume is a derived quantity.
- Achieve the target volume distribution by **super-sampling and down-selection**, not
  by rejection inside the generation loop.
- Make per-rank work **fixed and independent** for clean MPI scaling.

### 2.1 High-level workflow

```
0. Read molecule; detect molecular point group  ............... (§2.2, §2.3)
1. Determine compatible (layer group, Wyckoff position) buckets
   symbolically from Z + site symmetry vs molecular symmetry .. (§2.4)
   For each bucket: alignment rotations + residual DOF ........ (§2.5)
2. Compute area window [A_min, A_max] ......................... (§2.6.3)
   Enumerate distinct epitaxy sublattices (HNF), reduce,
   bucket by 2D lattice class ................................. (§2.6)
3. For each (layer group, Wyckoff) bucket:
   target N structures, generate I × N raw structures
   (split evenly across MPI ranks):
     a. draw a lattice-independent "structure seed":
        Wyckoff free coordinates + molecular orientation ...... (§2.5.3)
     b. draw an epitaxy matrix from the bucket's class list ... (§2.6)
     c. draw the allowed in-plane offset c_par of the c vector  (§2.7.1)
     d. compute the out-of-plane height h analytically from
        the interlayer contact condition; c = c_par + (0,0,h) . (§2.7.3)
     e. realize Cartesian/fractional coordinates; check
        intermolecular distances with the vdW cutoff matrix;
        on failure, redraw (b)-(d) up to max_retry times ...... (§2.8)
4. Down-select I × N -> N per bucket so the realized volumes
   follow N(V_pred, sigma); extensible to other criteria ...... (§2.9)
5. Collect on rank 0, write geometry.out with metadata ........ (§2.10)
```

### 2.2 Molecular point-group detection

A new standalone module (own source file; the workflow only calls it) detects the full
molecular point group from the geometry. Coordinates are first shifted to the centroid

$$\bar{\mathbf r} = \tfrac1N \sum_i \mathbf r_i, \qquad \mathbf r_i \leftarrow \mathbf r_i - \bar{\mathbf r},$$

which is invariant under every molecular symmetry operation (it permutes atoms of equal
element), so all symmetry elements pass through the origin afterwards.

#### 2.2.1 Atom equivalence classes

Two atoms $i, j$ are *equivalent* iff they have the same element and the same distance
environment: writing the multiset of pairs $\mathrm{sig}(i) = \{(\mathrm{elem}_k,\, d_{ik})\}_{k \ne i}$
sorted lexicographically, atoms are equivalent iff their signatures match entrywise
within the tolerance $\tau$. Every symmetry operation permutes atoms **within**
equivalence classes; this both restricts candidate operations and accelerates matching.

#### 2.2.2 Candidate axes and completeness

Candidate axes (used both as rotation axes and as mirror normals) are the union of:

- the three principal axes of the gyration tensor
  $\;T = \sum_i \mathbf r_i \mathbf r_i^{\mathsf T}$ (symmetric 3×3, solved by Jacobi
  iteration);
- per equivalence class, for atoms $i, j$ in the class:
  $\;\mathbf r_i$, $\;\mathbf r_i + \mathbf r_j$, $\;\mathbf r_i - \mathbf r_j$,
  $\;\mathbf r_i \times \mathbf r_j$ (normalized, deduplicated, sign-canonicalized).

*Completeness argument.* Any rotation axis of the molecule is a principal axis of every
invariant symmetric tensor, in particular of $T$; when $T$ is non-degenerate this already
yields all possible axes. In degenerate (symmetric-top / spherical-top) cases the
class-derived candidates cover the remaining possibilities: an axis either passes through
class atoms ($\mathbf r_i$), bisects pairs ($\mathbf r_i + \mathbf r_j$), or is normal to
the plane of a pair ($\mathbf r_i \times \mathbf r_j$); mirror normals are differences
$\mathbf r_i - \mathbf r_j$ of swapped atoms or principal axes. Group closure (§2.2.5)
is the final safety net: any operation expressible as a product of found operations is
recovered even if its axis was missed.

#### 2.2.3 Trial operations and verification

For each candidate axis $\hat{\mathbf a}$, the following orthogonal matrices are tested
($n$ up to the largest equivalence-class size, capped):

- proper rotations (Rodrigues form)
  $$ R(\hat{\mathbf a}, \theta) = \cos\theta\, I + \sin\theta\, [\hat{\mathbf a}]_\times + (1-\cos\theta)\, \hat{\mathbf a}\hat{\mathbf a}^{\mathsf T}, \qquad \theta = 2\pi/n,\; n = 2, 3, \dots $$
- mirrors $\;\sigma(\hat{\mathbf n}) = I - 2\,\hat{\mathbf n}\hat{\mathbf n}^{\mathsf T}$,
- improper rotations $\;S(\hat{\mathbf a}, \theta) = R(\hat{\mathbf a}, \theta)\,\sigma(\hat{\mathbf a})$ for even $n \ge 4$,
- and globally the inversion $-I$ and identity $I$.

A trial matrix $M$ is **valid** iff there exists a bijection (permutation) $\pi$ of atoms
with matching elements such that

$$ \max_i \; \lVert M\mathbf r_i - \mathbf r_{\pi(i)} \rVert \;\le\; \tau . $$

The matching uses nearest-neighbor assignment with a used-flag (sufficient because
$\tau \ll$ half the minimum interatomic distance). The achieved maximum deviation is
stored per operation as a quality diagnostic — it measures how far the actual (e.g.
DFT-relaxed) geometry deviates from exact symmetry. The tolerance $\tau$ is a user
parameter (default 0.1 Å).

#### 2.2.4 Exact operation identity: (permutation, determinant)

A central design decision: an operation is **identified by the pair $(\pi, \det M)$**,
not by its matrix. For a rigid, nonlinear molecule this pair is exactly unique:

- if two valid operations $g, h$ share $(\pi, \det)$, then $g h^{-1}$ fixes every atom
  and has determinant +1; a nontrivial proper rotation fixing all atoms requires all
  atoms collinear (excluded: nonlinear), hence $g = h$;
- the determinant is required because, for a planar molecule, the molecular-plane mirror
  induces the *identity permutation* (every atom fixed) — the same permutation as $E$ —
  and similarly an in-plane $C_2$ and the perpendicular mirror containing it induce equal
  permutations. The determinant separates these pairs.

Consequences: deduplication and group closure become **exact integer operations** (no
matrix tolerance), and the group multiplication table is computed by permutation
composition, $\pi_{gh} = \pi_g \circ \pi_h$, $\det_{gh} = \det_g \det_h$.

(Linear molecules are detected separately via the gyration spectrum — two vanishing
eigenvalues — and labeled $C_{\infty v}$ / $D_{\infty h}$.)

#### 2.2.5 Matrix refinement and group closure

For each valid operation, the matrix is **refined** from its permutation by the
orthogonal Procrustes problem: with covariance
$H = \sum_i \mathbf r_{\pi(i)} \mathbf r_i^{\mathsf T}$ and SVD $H = U \Sigma V^{\mathsf T}$,

$$ M^\star \;=\; U\, \mathrm{diag}\big(1,\, 1,\, s\big)\, V^{\mathsf T}, \qquad s = \pm 1 \text{ chosen so } \det M^\star = \det M, $$

i.e. the orthogonal matrix that best maps each atom onto its image in the least-squares
sense. For planar molecules $H$ is rank-2; the third singular direction is completed by
the cross product of the first two, with the sign fixed by the stored determinant. This
yields clean, near-machine-precision matrices from imperfect geometries.

The operation set is then closed under multiplication (in permutation space, exactly),
products refined as above, until no new $(\pi, \det)$ appears. The final count must
equal a valid point-group order; the per-operation deviations are reported.

#### 2.2.6 Point-group classification

From the inventory of elements, the Schoenflies symbol follows the standard flowchart:
linear → $C_{\infty v}$/$D_{\infty h}$ (by inversion); two or more axes of order ≥ 3 →
cubic family ($C_5$ → $I/I_h$; $C_4$ → $O/O_h$; else $T/T_d/T_h$); otherwise with
principal order $n$: $n$ perpendicular $C_2$ axes → $D_n / D_{nh} / D_{nd}$ (by
$\sigma_h$ / $n\,\sigma_d$); else $C_n / C_{nh} / C_{nv} / S_{2n}$; no proper axis →
$C_s / C_i / C_1$.

Expected results for the project's example molecules: naphthalene → $D_{2h}$ (order 8),
PTCDA → $D_{2h}$ (order 8); both planar within ~0.01 Å. The module is unit-tested on
these two geometries plus exact synthetic cases that exercise paths the planar $D_{2h}$
molecules do not: benzene ($D_{6h}$, degenerate top, order-6 axis), ammonia ($C_{3v}$,
axis not obtainable from atom pairs alone), methane ($T_d$, spherical top), and a chiral
molecule ($C_1$).

### 2.3 Representation of symmetry elements

Because the molecule is recentered, **every element passes through the origin**; elements
are therefore marked by directions and integer orders only — no positions needed:

| Element | Marking | Geometric meaning |
|---|---|---|
| identity $E$ | — | — |
| inversion $i$ | flag | point = centroid |
| proper axis $C_n$ | $(\hat{\mathbf a}, n)$, $\lVert\hat{\mathbf a}\rVert = 1$ | rotations by $2\pi k/n$ about the line $\{t\,\hat{\mathbf a}\}$ |
| mirror $\sigma$ | unit normal $\hat{\mathbf n}$ | plane $\{\mathbf r : \hat{\mathbf n}\cdot\mathbf r = 0\}$ |
| improper axis $S_n$ | $(\hat{\mathbf a}, n)$ | $R(\hat{\mathbf a}, 2\pi/n)$ followed by $\sigma(\hat{\mathbf a})$ |

So a 4-fold rotation axis is marked by its unit vector and the integer 4; a mirror plane
by its unit normal. Axes are sign-canonicalized ($\hat{\mathbf a}$ and $-\hat{\mathbf a}$
are the same element; the representative has its largest-magnitude component positive).
In addition to the element list, every individual operation is stored as its explicit
3×3 orthogonal Cartesian matrix (refined, §2.2.5). Useful extras recorded by the module:
planarity flag and molecular-plane normal, principal axes and moments, per-operation
maximum deviation.

### 2.4 Compatible (layer group, Wyckoff position) determination — symbolic

#### 2.4.1 Compatibility condition

No fake lattice and no numeric overlap test. A Wyckoff position of a layer group is
compatible iff:

1. **Multiplicity**: the position's multiplicity equals Z; and
2. **Oriented subgroup condition**: there exists a rotation $R \in SO(3)$ with

$$ R\, M\, R^{-1} \;\supseteq\; S, $$

where $M$ is the molecular point group (as Cartesian matrices, §2.2) and $S$ is the
site-symmetry group of the position (as Cartesian matrices in the crystal frame, §2.5.1).
This is a *conjugacy* condition, not abstract isomorphism: improper site operations
(mirrors, $\bar n$) can only be matched by improper molecular operations — a chiral
molecule cannot occupy a mirror site. Since crystallographic site groups are small, the
check is constructive (§2.5.2) and exact.

#### 2.4.2 Stabilizer / overlap list, computed symbolically

Let $Z_{\mathrm{gen}}$ denote the multiplicity of the general position = the number of
layer-group operations. When a molecule occupies a special position of multiplicity Z and
all $Z_{\mathrm{gen}}$ operations are applied, each distinct molecule appears
$Z_{\mathrm{gen}}/Z$ times. The **overlap list** — which generated copies coincide with
the first molecule — is the site stabilizer, computed exactly: operation $(W, \mathbf w)$
stabilizes the representative coordinate $\mathbf x_0$ iff

$$ W \mathbf x_0 + \mathbf w \;\equiv\; \mathbf x_0 \pmod{L_{2D}} \quad (\text{z compared exactly}), $$

with $L_{2D}$ the in-plane lattice translations. The coset decomposition of the group by
this stabilizer gives the mapping from operations to distinct molecules, replacing the
current on-the-fly numeric overlap detection. It is used to verify placements, merge
coincident copies, and tell the distance checker which image pairs are "self".

#### 2.4.3 Buckets, deduplication, naming

- Multiple Wyckoff positions per layer group can be compatible simultaneously (e.g.
  layer group p$\bar 1$ has four distinct inversion-center positions of multiplicity 1
  plus the general position of multiplicity 2). The generation unit is therefore the
  **(layer group, Wyckoff position) bucket**, not the layer group.
- Positions equivalent under the affine normalizer of the group (e.g. the four inversion
  centers of p$\bar 1$, differing only by origin shift) generate identical films up to a
  translation that is later absorbed by the film–substrate registry; one representative
  per normalizer orbit is kept.
- The per-bucket quota parameter is named `num_structures_per_wyckoff_position`
  (replacing the per-layer-group `num_structures`). Total output =
  $N \times \#\text{buckets}$; the bucket list is printed at startup.

### 2.5 Molecule placement on a Wyckoff position

#### 2.5.1 The ideal crystal frame (lattice-independent)

Define, per lattice class, the orthonormal frame
$$ \hat{\mathbf z} = \text{out-of-plane unit vector}, \qquad \hat{\mathbf x} = \mathbf a / \lVert \mathbf a \rVert, \qquad \hat{\mathbf y} = \hat{\mathbf z} \times \hat{\mathbf x}. $$

In this frame the site-symmetry operations of every Wyckoff position are **fixed
matrices**, independent of the lattice parameters: oblique classes reference only
$\hat{\mathbf z}$; rectangular classes have 2-folds and mirror normals along
$\hat{\mathbf x}, \hat{\mathbf y}$; square adds the diagonals; hexagonal uses
$\hat{\mathbf x}$ rotated by multiples of 30°/60°. This is what makes molecular
orientation representable *before* any lattice exists.

#### 2.5.2 Alignment rotations

For a compatible bucket, the alignment rotations are computed once:

1. Choose generators of the site group $S$ (at most two suffice for crystallographic
   site groups).
2. For each generator (kind, order, axis $\hat{\mathbf u}$ in the ideal frame), select a
   molecular element of the same kind and order with axis $\hat{\mathbf m}$; a first
   match constrains $R\,\hat{\mathbf m} = \pm\hat{\mathbf u}$, determining $R$ up to a
   spin about $\hat{\mathbf u}$; a second independent match fixes $R$ completely.
3. Verify the full condition $R M R^{-1} \supseteq S$ operation by operation.
4. Collect all distinct solutions $R_0^{(1)}, R_0^{(2)}, \dots$ **modulo proper rotations
   of $M$ itself** (two alignments differing by a proper molecular symmetry produce the
   same physical placement). The result is a finite list of discrete orientation choices.

#### 2.5.3 Residual degrees of freedom and the structure seed

After alignment, the remaining *continuous* orientation freedom is the centralizer of
$S$ in $SO(3)$ — rotations commuting with every site operation:

| Site symmetry | Continuous orientational DOF | Description |
|---|---|---|
| $1$, $\bar 1$ | 3 | fully free orientation (random uniform rotation) |
| $2$, $m$, $2/m$, $3$, $4$, $6$, $\bar 3$, $\bar 4$, $\bar 6$, $4/m$, $6/m$ | 1 | spin angle $\varphi$ about the matched axis |
| $222$, $mm2$, $mmm$, $3m$, $4mm$, $\bar42m$, $6mm$, … | 0 | only the discrete alignment choices of §2.5.2 |

This reproduces and generalizes the existing 2/1/0 classification
(`lg_get_degrees_of_freedom`) — a regression anchor for the new code. Highly restricted
positions (0 DOF) are acceptable: the bucket still produces $I \times N$ structures,
with diversity supplied by the epitaxy matrix, the c-vector offset, and the Wyckoff free
coordinates.

A raw structure is represented by the lattice-independent **structure seed**:

$$ \big(\;\text{layer group},\ \text{Wyckoff position},\ \text{alignment choice } k,\ x_{\mathrm{frac}},\, y_{\mathrm{frac}},\ z_{\mathrm{cart}},\ \varphi \text{ or quaternion}\;\big) $$

- $x_{\mathrm{frac}}, y_{\mathrm{frac}}$: the position's free in-plane coordinates,
  sampled uniformly in $[0,1)$ (fixed coordinates taken from the Wyckoff representative).
- $z_{\mathrm{cart}}$: the COM height in **Cartesian** Å (not fractional — fractional z
  would depend on the not-yet-determined cell height). For site symmetries containing
  any z-flipping element, $z$ is pinned to 0 (in layer groups all such elements lie in
  the $z=0$ plane); otherwise it parameterizes the molecule's height within the film and
  is sampled within the molecule's extent.
- Atomic *fractional* coordinates are computed only at realization time (§2.8), because
  the molecule is a rigid Cartesian body: its fractional coordinates depend on the
  lattice, but its COM-and-orientation description does not.

At realization, the actual film vectors give a single in-plane rotation $T$ mapping the
ideal frame onto the actual cell directions, and the Cartesian orientation is
$R_{\mathrm{cart}} = T \cdot R_{\mathrm{residual}}(\varphi) \cdot R_0^{(k)}$.

### 2.6 Epitaxy matrices

#### 2.6.1 Canonical enumeration (Hermite Normal Form)

The film's in-plane lattice must be a sublattice of the substrate surface lattice
$(\mathbf s_1, \mathbf s_2)$. Distinct sublattices of index $d$ are in one-to-one
correspondence with integer **Hermite Normal Form** matrices

$$ H = \begin{pmatrix} p & q \\ 0 & r \end{pmatrix}, \qquad p\,r = d, \quad 0 \le q < r \;\; (p, r \ge 1), $$

so enumerating HNFs enumerates each sublattice **exactly once** — eliminating the
unimodular redundancy of the current brute-force scan (any
$U H$ with $U \in GL(2,\mathbb Z)$, $\det U = 1$, describes the same film lattice). The
number of sublattices of index $d$ is the divisor sum $\sigma(d)$; summed over the area
window below this is typically $10^2$–$10^4$ candidates — kilobytes instead of the
current 1.6 GB, enumerated in milliseconds.

Random on-the-fly generation of integer matrices was considered and rejected: it samples
sublattices with multiplicity proportional to their number of representations (a silent
bias), produces duplicate films, cannot detect infeasible layer groups up front, and has
no cost advantage once the area window bounds the candidate set.

#### 2.6.2 Reduction and lattice-class bucketing

Each HNF basis is Lagrange (Gauss) reduced to the shortest in-plane basis; from it,
$(\lVert\mathbf a\rVert, \lVert\mathbf b\rVert, \gamma, \text{area})$ are computed and
the candidate is bucketed by 2D lattice class with tolerances (γ = 90° → rectangular;
$\lVert\mathbf a\rVert = \lVert\mathbf b\rVert$, γ = 90° → square; γ = 120° equal
lengths → hexagonal; else oblique), subject to the constraints: right-handedness
(det > 0), aspect ratio $\lVert\mathbf a\rVert / \lVert\mathbf b\rVert < 3$, angle
$\gamma \in [30°, 150°]$, and the area window of §2.6.3. The epitaxial constraint means
classes are only realized as exactly as the substrate permits; the bucketing tolerance is
an explicit parameter, as in the current code. At generation time a matrix is sampled
uniformly (optionally area-stratified) from the bucket matching the layer group's class.
Buckets that are empty for a given class mark the corresponding layer groups infeasible
**before** any generation attempt.

#### 2.6.3 Area window

- **Hard floor.** From the molecule's principal-axes-aligned vdW bounding box: along
  each principal axis $\hat{\mathbf e}_k$,
  $\;\mathrm{ext}_k = \max_i(\mathbf r_i\cdot\hat{\mathbf e}_k + r_{\mathrm{vdw},i}) - \min_i(\mathbf r_i\cdot\hat{\mathbf e}_k - r_{\mathrm{vdw},i})$;
  the minimal facet is the product of the two smallest extents, and
  $A_{\min} = \alpha \cdot \mathrm{ext}_{(1)}\mathrm{ext}_{(2)}$ with safety factor
  $\alpha \approx 0.7$ (the rectangle overestimates the true silhouette, and non-convex
  molecules can interdigitate). Note the floor uses **one** molecule, not Z of them:
  layer-group operations (mirrors, inversion, in-plane 2-folds) may place the Z molecules
  at different heights, so their in-plane footprints can legally overlap;
  $Z\times$cross-section would wrongly exclude stacked packings. The $Z$-scaled value is
  retained only as a soft diagnostic.
- **Ceiling.** $A_{\max} = V_{\mathrm{pred}} / h_{\mathrm{floor}}$ with
  $h_{\mathrm{floor}} \approx 2$ Å as the minimal film height (equivalently the current
  $0.5 \times V$ rule), used as a fallback bound.
- **Volume pre-windowing (efficiency, not correctness).** Since volume is derived as
  $V = A \cdot h$ (§2.7.4) and the film thickness is bracketed by the molecular extents,
  $h \in [t_{\min}, t_{\max}]$, restricting
  $\;A \in [(V_{\mathrm{pred}} - k\sigma)/t_{\max},\ (V_{\mathrm{pred}} + k\sigma)/t_{\min}]$
  pre-aligns the raw pool with the target volume Gaussian so that the down-selection of
  §2.9 discards little.

### 2.7 The out-of-plane vector c

#### 2.7.1 Direction: derived from symmetry closure, not hardcoded

Write $\mathbf c = \mathbf c_\parallel + (0, 0, h)$ with $\mathbf c_\parallel$ the
in-plane component. The layer-group operations remain symmetry operations of the
3D-stacked crystal iff for every operation with point part $W$:

$$ W \mathbf c \;\in\; \varepsilon_W\, \mathbf c + L_{2D}, \qquad \varepsilon_W = \begin{cases} +1 & W \text{ preserves } z \\ -1 & W \text{ flips } z \end{cases} $$

Evaluating this rule per layer group over its operation list yields the allowed
$\mathbf c_\parallel$ (this generalizes the per-class treatment of the current code and
reproduces its continuous tilts):

| Layer groups | Constraint on $\mathbf c_\parallel$ |
|---|---|
| 1–2 (p1, p$\bar 1$) | unconstrained — 2 continuous parameters |
| 3–7 (z-axis operations) | $2\mathbf c_\parallel \in L_{2D}$ → discrete: $\{0, \tfrac{\mathbf a}2, \tfrac{\mathbf b}2, \tfrac{\mathbf a + \mathbf b}2\}$ |
| 8–18 (unique in-plane axis ∥ x) | x-component discrete $\{0, \tfrac a2\}$; y-component continuous |
| 19–48, square, hexagonal | both components discrete; the allowed subset per group follows from the rule |

The discrete half-lattice offsets are the staggered (AB-type) stackings — symmetry-legal
and physically important for the planned **multilayer generation by replication along
c**; the current code omits them (it uses only $\mathbf c_\parallel = 0$ outside the two
continuous classes). The new scheme samples $\mathbf c_\parallel$ from the allowed set
(uniform over discrete options; uniform over the reduced cell for continuous components)
as one of the pool's degrees of freedom.

#### 2.7.2 Implementation: generate perpendicular, then tilt

A tilted c changes the *fractional* matrices of the symmetry operations (their third
column acquires integer entries), which would ripple through all symmetry-application
code. Instead: **generate the film with c ⊥ the ab-plane** (database operations apply
unchanged), then set the final cell vector $\mathbf c = \mathbf c_\parallel + (0,0,h)$
while keeping the atoms' Cartesian positions, and re-express fractional coordinates once
at output. The 2D layer is untouched — only the stacking registry of its periodic images
changes — and the closure rule of §2.7.1 guarantees the layer-group operations remain
crystal symmetries of the result.

#### 2.7.3 Height: analytic interlayer contact condition

The height is computed **analytically** — no rejection loop — from the requirement that
the closest interlayer atomic contact equals a prescribed multiple of the vdW radii sum.
Let atoms have Cartesian positions $\mathbf r_i = (x_i, y_i, z_i)$ (all Z molecules,
after realization on the chosen in-plane lattice), and prescribe per-pair minimum
distances $D_{ij} = f_z\,(r_{\mathrm{vdw},i} + r_{\mathrm{vdw},j})$ with $f_z$ a user
parameter (default 1.5). For the image shifted by
$\mathbf c = \mathbf c_\parallel + (0,0,h)$ plus any in-plane translation
$\mathbf t \in L_{2D}$, the constraint
$\lVert \mathbf r_i - \mathbf r_j - \mathbf t - \mathbf c \rVert \ge D_{ij}$ depends on
$h$ only through the z-gap. With the in-plane image separation
$\rho_{ij,\mathbf t} = \lVert (\Delta x_{ij}, \Delta y_{ij}) - \mathbf t_{xy} - \mathbf c_\parallel \rVert$:

- pairs with $\rho \ge D_{ij}$ impose nothing;
- pairs with $\rho < D_{ij}$ require $\;|z_i - z_j - h| \ge \sqrt{D_{ij}^2 - \rho^2}$,

giving the exact minimal height in closed form:

$$ \boxed{\; h_{\min} \;=\; \max_{\substack{i,\,j,\ \mathbf t \,:\, \rho_{ij,\mathbf t} < D_{ij}}} \Big( z_i - z_j + \sqrt{D_{ij}^2 - \rho_{ij,\mathbf t}^2} \Big) \;} $$

The maximum runs over $O((NZ)^2)$ pairs times the few in-plane images with
$\rho < \max D_{ij}$ (the image search radius is a few Å). Setting $h = h_{\min}$ puts
the limiting interlayer pair exactly at distance $D_{ij}$ — true 3D distance, not just
z-gap — and every other pair at or beyond its own threshold. Intra-layer distances are
independent of $h$ and are handled by the in-plane check (§2.8). The formula is the
rigid-body "drop the layer onto its own periodic image" contact problem; it already
accommodates the $\mathbf c_\parallel$ offset of §2.7.1, so the direction decision and
the height calculation are independent.

#### 2.7.4 Why z-first, volume-derived (design decision)

The height is **not** chosen to meet a target volume. Rationale: films will later be
replicated along c for multilayer generation; padding h beyond contact inserts vacuum
between layers, producing structures that are not meaningful packings and that break
multilayer construction. With $h = h_{\min}$, every structure in the pool is a physically
sensible stacking, the interlayer spacing is controlled explicitly by $f_z$, and the
cell volume is a *derived* quantity

$$ V = A_{\mathrm{cell}} \cdot h_{\min}, $$

(the in-plane offset $\mathbf c_\parallel$ does not change the volume). The target
volume distribution is then imposed by down-selection (§2.9), assisted by the area
pre-windowing of §2.6.3.

### 2.8 Realization and acceptance check

For each raw structure: realize Cartesian coordinates (in-plane lattice from the epitaxy
matrix; orientation via the frame rotation of §2.5.3; symmetry copies via the layer-group
operations with perpendicular c; stabilizer-based merging per §2.4.2), compute
$h_{\min}$, assemble the final cell, and check all intermolecular atomic distances —
including in-plane periodic images — against the vdW cutoff matrix (as today). By
construction the $\pm\mathbf c$ images cannot fail (the contact condition uses
$f_z \ge$ the check's factor), so the check is effectively in-plane. On failure, redraw
the epitaxy matrix and $\mathbf c_\parallel$ (keeping the structure seed) up to
`max_retry` times (user parameter, default 5000); a seed that exhausts its retries is
regenerated. Per-bucket attempt and acceptance counts are recorded.

### 2.9 Super-sampling and down-selection

With user quota $N$ = `num_structures_per_wyckoff_position` and integer super-sampling
ratio $I \ge 1$ (`super_sampling_ratio`, default 10), each bucket generates $I \times N$
accepted raw structures, then down-selects to $N$ so the realized volumes follow the
target Gaussian $\mathcal N(V_{\mathrm{pred}}, \sigma^2)$. Selection is by importance
resampling: weight structure $s$ by

$$ w_s \;\propto\; \frac{\mathcal N(V_s;\, V_{\mathrm{pred}}, \sigma^2)}{\hat p(V_s)}, $$

with $\hat p$ the empirical (e.g. histogram/KDE) density of the raw pool, then sample $N$
without replacement. This shapes the distribution only where the pool has support — the
area pre-windowing (§2.6.3) is what ensures adequate coverage. The framework extends to
future criteria (packing compactness, stacking quality, diversity filters) as additional
multiplicative weights.

### 2.10 MPI parallelization

Work is fixed and embarrassingly parallel: each of the $P$ ranks generates
$I \times N / P$ raw structures per bucket with rank-decorrelated seeds, with **no
synchronization during generation** (no per-batch gather/broadcast, no early-stop
coordination). Ranks send accepted raw structures (or, cheaper, their seeds + epitaxy
choice + derived scalars) to rank 0 at bucket completion; rank 0 performs the
down-selection of §2.9 and writes `geometry.out` with the metadata required downstream:
layer group, Wyckoff position, epitaxy matrix, COM, orientation (Euler angles /
quaternion), $\mathbf c_\parallel$ choice, realized volume, contact pair info.

### 2.11 Parameters

| Parameter | Status | Note |
|---|---|---|
| `vdw_matrix` | kept | acceptance check |
| `Z` | kept | |
| `volume_mean`, `volume_std` | kept | target Gaussian for down-selection ($V_{\mathrm{pred}}, \sigma$) |
| `lattice_vector_2d_from_geo` | kept | substrate surface lattice |
| `random_seed` | kept | |
| `tol` | kept (re-scoped) | now the symmetry-detection tolerance $\tau$ |
| `num_structures` | renamed | → `num_structures_per_wyckoff_position` ($N$, per bucket) |
| `max_attempts` | repurposed | → `max_retry` per structure (default 5000) |
| `interface_area_mean/std` | deprecated | replaced by derived area window (§2.6.3); optional override |
| `volume_multiplier` | deprecated | replaced by analytic $h$ |
| `vol_attempt` | deprecated | no volume-resampling cadence |
| `spg_dist_type` | deprecated | uniform-per-bucket implied; multiplicity weighting can be reintroduced on $N$ if needed |
| `Zp_max`, `sr` | deprecated | already unused in the layer path |
| `super_sampling_ratio` | **new** | $I \ge 1$, integer, default 10 |
| `z_contact_factor` | **new** | $f_z$, default 1.5 |
| internal knobs | new | area-floor safety $\alpha$, lattice-class bucketing tolerance, $k\sigma$ volume window |

### 2.12 Validation and regression plan

1. **Symmetry module**: unit tests on naphthalene and PTCDA (both $D_{2h}$, order 8,
   planar) and on exact synthetic molecules (benzene $D_{6h}$, ammonia $C_{3v}$, methane
   $T_d$, chiral $C_1$); per-operation deviations reported.
2. **Compatibility module**: diff the symbolic (layer group, Wyckoff) list and DOF
   classification against the current fake-lattice detector and
   `lg_get_degrees_of_freedom` on the example molecules.
3. **Epitaxy module**: verify HNF counts against $\sigma(d)$; verify reduced bases
   reproduce the current sampler's accepted combos; confirm class buckets.
4. **c-vector / height**: brute-force numeric check of $h_{\min}$ (scan h, verify
   contact pair and feasibility boundary); closure-rule results spot-checked against the
   current per-class c constructions.
5. **End-to-end**: generated structures re-checked with the existing
   vdW-matrix checker; symmetry of outputs verified (layer group detection); volume
   histogram vs target Gaussian.

### 2.13 Module decomposition and implementation order

1. `molecule_symmetry` — point-group detection (§2.2–2.3), standalone file + tests.
2. `site_compatibility` — symbolic bucket determination, alignment rotations, stabilizer
   lists (§2.4–2.5), validated against the current detector.
3. `epitaxy_enumeration` — HNF enumeration, reduction, class bucketing, area window
   (§2.6).
4. `c_vector` — closure rule for $\mathbf c_\parallel$, analytic $h_{\min}$ (§2.7).
5. New generation driver — seeds, realization, acceptance, super-sampling pool,
   down-selection, MPI collection, output (§2.8–2.10), replacing the body of
   `mpi_generate_layer_with_vdw_cutoff_matrix` while preserving its output format and
   downstream metadata.

Each stage has a test oracle before the next depends on it.

### 2.14 Open questions

- Exact default of the area-floor safety factor $\alpha$ and the volume pre-window
  width $k$ (tune on the example systems).
- Whether to expose area-window overrides to the user or keep them fully derived.
- Down-selection extensions (compactness, stacking-quality metrics) — deferred by design.
- Multilayer generation itself (replication count, inter-layer registry optimization) —
  future work building on §2.7.
