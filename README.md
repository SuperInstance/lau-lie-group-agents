# lau-lie-group-agents

**Lie group theory for autonomous agents** — SO(n), SU(n), GL(n), Sp(2n), SE(3), Lie algebras with the bracket, the exponential map, Baker-Campbell-Hausdorff formula, adjoint representation, structure constants, Killing form, root systems (Aₙ, Bₙ, Cₙ, Dₙ), and the Peter-Weyl theorem.

## What This Does

This crate provides concrete, computational implementations of the foundational structures of Lie theory — the mathematical language of continuous symmetry. Every structure is represented as a real matrix using `nalgebra`, making it immediately usable for robotics, physics simulations, and agent-based systems that need to reason about rotations, rigid body motions, and symmetry groups.

You get:

- **Five classical Lie groups**: SO(n), SU(n), GL(n), Sp(2n), SE(3) — with membership tests, constructors, multiplication, and inverses
- **Lie algebra operations**: the bracket [X, Y] = XY − YX, basis construction for so(n) and sp(2n), Jacobi identity verification
- **Exponential map**: matrix exponential (scaling-and-squaring + Taylor series) and matrix logarithm, mapping algebra → group and back
- **BCH formula**: log(exp(X)exp(Y)) at first, second, and third order
- **Adjoint representation**: Ad(g) and ad(X) as matrices, with Ad(exp X) = exp(ad X) verification
- **Structure constants**: the tensor c^k_{ij} for so(n), with antisymmetry and Jacobi checks
- **Killing form**: B(X, Y) = Tr(ad X ∘ ad Y), with negative-definiteness, bilinearity, symmetry, and invariance tests
- **Root systems**: Aₙ, Bₙ, Cₙ, Dₙ with positive roots, simple roots, and Cartan matrices
- **Peter-Weyl theorem**: representation decomposition, characters, and the Weyl character formula for SU(2)

## Key Idea

A **Lie group** is a group that is also a smooth manifold — you can multiply elements and take inverses smoothly. The tangent space at the identity is the **Lie algebra**, a vector space with an antisymmetric bilinear operation (the bracket) satisfying the Jacobi identity. The exponential map bridges the two: exp: g → G sends algebra elements to group elements.

For agents, Lie groups describe the space of possible transformations: rotations (SO(3)), rigid body motions (SE(3)), and more exotic symmetries. The Lie algebra is the space of "infinitesimal" transformations — velocities, angular momenta, forces. The bracket measures the failure of transformations to commute, which is the fundamental obstruction to integrating agent behaviors in parallel.

## Install

```toml
[dependencies]
lau-lie-group-agents = "0.1.0"
```

Requires `nalgebra` (with `serde-serialize`), `serde` + `serde_json`, and `num-complex`.

## Quick Start

```rust
use lau_lie_group_agents::lie_groups::{SO, SE3, Sp};
use lau_lie_group_agents::lie_algebra::{SoAlgebra, LieAlgebraElement, verify_jacobi};
use lau_lie_group_agents::exponential_map::matrix_exp;
use lau_lie_group_agents::bch::{bch, BCHOrder};

// SO(3) rotation around z-axis by 90°
let r = SO::rotz(std::f64::consts::FRAC_PI_2);
assert!(SO::is_member(&r.matrix));

// Lie algebra basis for so(3): 3 antisymmetric matrices
let basis = SoAlgebra::basis(3);
assert_eq!(basis.len(), 3); // dim so(3) = 3

// Bracket [e₁, e₂] gives another basis element
let e01 = &basis[0];
let e12 = &basis[1];
let bracket = LieAlgebraElement::bracket(e01, e12);
assert!(bracket.is_antisymmetric());

// Exponential map: algebra → group
let rot = matrix_exp(&e01.matrix);
assert!(SO::is_member(&rot));

// BCH: log(exp(X)exp(Y)) ≈ X + Y + ½[X,Y] + ...
let bch_result = bch(e01, e12, BCHOrder::Third);

// Rigid body transform: rotate then translate
let t = SE3::translation(1.0, 2.0, 3.0);
let (rot, trans) = SE3::decompose(&t.matrix);
```

## API Reference

### `lie_groups` — Classical Matrix Lie Groups

| Type | Description |
|------|-------------|
| `LieGroupElement` | A group element stored as a matrix with a `GroupType` tag. |
| `GroupType` | Enum: `SO(n)`, `SU(n)`, `GL(n)`, `Sp(n)`, `SE3`. |

**`LieGroupElement` methods:**

- `new(matrix, group_type)` — Create from a matrix.
- `identity(group_type)` — The identity element.
- `multiply(&other)` — Group multiplication (matrix product).
- `inverse() → Option<Self>` — Matrix inverse (fails if singular).
- `dim() → usize` — Matrix dimension.

**`GroupType::matrix_dim() → usize`** — Returns the matrix dimension.

#### SO(n) — Special Orthogonal Group

Accessed via `SO::` static methods.

- `SO::identity(n)` — Identity matrix.
- `SO::rotation(n, i, j, theta)` — Rotation by `theta` in the (i,j)-plane.
- `SO::is_member(&m) → bool` — Checks M^T M = I and det M = 1.
- `SO::so2(theta)` — 2D rotation.
- `SO::rotx(theta)`, `SO::roty(theta)`, `SO::rotz(theta)` — SO(3) axis rotations.

#### GL(n) — General Linear Group

- `GL::identity(n)`, `GL::from_matrix(n, m)`, `GL::diagonal(n, &[entries])`
- `GL::is_member(&m) → bool` — Checks det ≠ 0.

#### SU(n) — Special Unitary Group (realified)

- `SU::identity(n)`, `SU::su2_rotation(theta)` — SU(2) element as 2×2 real matrix.
- `SU::is_member_real(&m) → bool` — Checks the realified unitary conditions.

#### Sp(2n) — Symplectic Group

- `Sp::identity(n)`, `Sp::j_matrix(n)` — The standard symplectic form J.
- `Sp::is_member(&m, n) → bool` — Checks M^T J M = J.
- `Sp::symplectic_rotation(n, theta)` — Block-diagonal symplectic rotation.

#### SE(3) — Special Euclidean Group

- `SE3::identity()`, `SE3::from_rotation_translation(&rot, &trans)`
- `SE3::translation(tx, ty, tz)` — Pure translation.
- `SE3::decompose(&m) → (DMatrix, DVector)` — Extract rotation and translation.
- `SE3::transform_point(&m, &p) → DVector` — Apply rigid transform to a 3D point.

---

### `lie_algebra` — Lie Algebra Elements and Operations

| Type | Description |
|------|-------------|
| `LieAlgebraElement` | An n×n matrix representing an algebra element. |
| `SoAlgebra` | Factory for so(n) basis elements. |
| `SpAlgebra` | Membership test for sp(2n). |

**`LieAlgebraElement` methods:**

- `new(matrix)`, `zero(n)` — Constructors.
- `bracket(&x, &y)` — The Lie bracket [X, Y] = XY − YX.
- `add(&other)`, `scale(s)` — Vector space operations.
- `dim() → usize`, `norm() → f64`
- `is_traceless() → bool` — For sl(n), su(n).
- `is_antisymmetric() → bool` — For so(n).

**`SoAlgebra`:**

- `basis_element(n, i, j)` — The basis element E_{ij} − E_{ji}.
- `basis(n) → Vec<LieAlgebraElement>` — Full basis (n(n−1)/2 elements).
- `dim(n) → usize` — Returns n(n−1)/2.

**`SpAlgebra`:**

- `is_member(&m, n) → bool` — Checks XJ + JX^T = 0.

**Free functions:**

- `verify_jacobi(&x, &y, &z) → bool` — Checks [X,[Y,Z]] + [Y,[Z,X]] + [Z,[X,Y]] = 0.
- `verify_antisymmetry(&x, &y) → bool` — Checks [X,Y] = −[Y,X].
- `verify_bilinearity(&x, &y, &z, a, b) → bool` — Checks [aX+bY, Z] = a[X,Z] + b[Y,Z].

---

### `exponential_map` — Matrix Exponential and Logarithm

| Type | Description |
|------|-------------|
| `ExponentialMap` | Configurable exponential map for a given dimension. |

**Free functions:**

- `matrix_exp(&m) → DMatrix<f64>` — Compute exp(M) via scaling-and-squaring with Taylor series.
- `matrix_log(&m) → DMatrix<f64>` — Compute log(M) via inverse scaling-and-squaring.

**`ExponentialMap` methods:**

- `new(n)` — Create for n×n matrices.
- `apply(&x) → LieGroupElement` — exp: g → G.
- `verify_identity() → bool` — exp(0) = I.
- `verify_so_mapping() → bool` — Checks exp maps so(n) → SO(n).
- `exp_t(&x, t) → DMatrix<f64>` — Compute exp(tX).

---

### `bch` — Baker-Campbell-Hausdorff Formula

| Type | Description |
|------|-------------|
| `BCHOrder` | Enum: `First` (X+Y), `Second` (+½[X,Y]), `Third` (+1/12 terms). |

**Free functions:**

- `bch(&x, &y, order) → LieAlgebraElement` — Compute log(exp(X)exp(Y)) to given order.
- `verify_bch(&x, &y, order) → bool` — Check exp(X)exp(Y) ≈ exp(BCH(X,Y)).
- `bch_series(&x, &y, depth) → LieAlgebraElement` — Compute to arbitrary depth.

---

### `adjoint` — Adjoint Representation

| Type | Description |
|------|-------------|
| `AdjointRep` | The adjoint representation Ad: G → GL(g) and ad: g → End(g). |

**Methods:**

- `new(n)` — Create for n-dimensional group.
- `ad_g(&g, &x) → LieAlgebraElement` — Ad(g)X = gXg⁻¹.
- `ad_x(&x, &y) → LieAlgebraElement` — ad(X)Y = [X,Y].
- `ad_matrix(&x, &basis) → DMatrix<f64>` — Matrix of ad(X) in the given basis.
- `verify_ad_exp(&x) → bool` — Check Ad(exp X) = exp(ad X).
- `trace_ad(&x, &basis) → f64` — Trace of ad(X) (0 for semisimple algebras).

---

### `structure_constants` — Structure Constants

| Type | Description |
|------|-------------|
| `StructureConstants` | The tensor c^k_{ij} defined by [eᵢ, eⱼ] = Σ c^k_{ij} eₖ. |

**Methods:**

- `from_basis(&[elements]) → Self` — Compute from any basis.
- `for_so(n) → Self` — Compute for so(n).
- `get(i, j, k) → f64` — Get c^k_{ij}.
- `verify_antisymmetry() → bool` — c^k_{ij} = −c^k_{ji}.
- `verify_jacobi() → bool` — Checks the Jacobi identity on structure constants.
- `dim() → usize`

---

### `killing_form` — Killing Form

| Type | Description |
|------|-------------|
| `KillingForm` | B(X, Y) = Tr(ad X ∘ ad Y), the natural invariant bilinear form. |

**Methods:**

- `for_so(n) → Self` — Create for so(n).
- `evaluate(&x, &y, &basis) → f64` — B(X, Y).
- `matrix(&basis) → DMatrix<f64>` — The full matrix B_{ij}.
- `is_negative_definite(n) → bool` — True for compact semisimple algebras.
- `verify_bilinearity(&basis) → bool`
- `verify_symmetry(&basis) → bool` — B(X,Y) = B(Y,X).
- `verify_invariance(&basis) → bool` — B([Z,X],Y) + B(X,[Z,Y]) = 0.

---

### `root_systems` — Root Systems Aₙ, Bₙ, Cₙ, Dₙ

| Type | Description |
|------|-------------|
| `RootSystemType` | Enum: `A`, `B`, `C`, `D`. |
| `RootSystem` | A complete root system with positive roots, simple roots, and Cartan matrix. |

**`RootSystem` static constructors:**

- `type_a(n)` — For su(n+1). Positive roots: eᵢ − eⱼ (i < j).
- `type_b(n)` — For so(2n+1). Adds short roots eᵢ.
- `type_c(n)` — For sp(2n). Adds long roots 2eᵢ.
- `type_d(n)` — For so(2n). eᵢ ± eⱼ only.

**Methods:**

- `num_positive_roots() → usize`
- `num_simple_roots() → usize`
- `num_roots() → usize` — Total (positive + negative).
- `cartan_matrix() → DMatrix<f64>` — C_{ij} = 2(αᵢ, αⱼ)/(αᵢ, αᵢ).
- `verify_cartan_diagonal() → bool` — Checks 2's on the diagonal.

**Free function:**

- `count_positive_roots(system_type, rank) → usize` — Formula: Aₙ→n(n+1)/2, Bₙ→n², Cₙ→n², Dₙ→n(n−1).

---

### `peter_weyl` — Peter-Weyl Theorem and Representation Theory

| Type | Description |
|------|-------------|
| `Representation` | A unitary representation with dimension and type. |
| `RepType` | Enum: `Trivial`, `Standard`, `Adjoint`, `Fundamental`, `Custom(n)`. |
| `PeterWeylDecomposition` | L²(G) decomposed into irreducible representations. |

**`Representation` methods:**

- `trivial()`, `standard(n)`, `adjoint(d)`, `fundamental(n)`
- `character(&g) → f64` — χ(g) = Tr(ρ(g)).
- `is_irreducible_heuristic() → bool` — Quick heuristic check.
- `dimension() → usize`

**`PeterWeylDecomposition` methods:**

- `for_so(n) → Self` — Standard decomposition: trivial + standard + adjoint.
- `num_representations() → usize`
- `total_dimension() → usize` — Σ dim(π)².
- `verify_orthogonality() → bool` — Simplified orthogonality check.
- `matrix_coefficient_inner_product(...)` — Peter-Weyl inner product formula.
- `regular_rep_dimension() → usize`

**Free functions:**

- `su2_character(n, theta) → f64` — Weyl character formula: χ_n(θ) = sin((n+1)θ)/sin(θ).
- `su_n_dimension(n, &highest_weight) → usize` — Weyl dimension formula.

## How It Works

All structures are represented concretely as real matrices:

1. **Group elements** are stored as their matrix representation (n×n for SO(n), 2n×2n for Sp(2n), 4×4 for SE(3)).

2. **Algebra elements** are also matrices — the Lie bracket is the matrix commutator [X, Y] = XY − YX.

3. **The exponential map** uses scaling-and-squaring: scale the matrix so ‖M/2^k‖ < 1, compute the Taylor series exp(M/2^k) ≈ I + M/2^k + (M/2^k)²/2! + …, then square k times.

4. **The logarithm** uses inverse scaling-and-squaring: repeatedly take element-wise square roots until close to identity, then apply the log Taylor series log(I + X) ≈ X − X²/2 + X³/3 − …

5. **Structure constants** are computed by expressing brackets in a chosen basis using the Frobenius inner product.

6. **Root systems** are constructed explicitly: for Aₙ, the positive roots are eᵢ − eⱼ (i < j); the Cartan matrix is computed from inner products of simple roots.

7. **The Killing form** B(X,Y) = Tr(ad X ∘ ad Y) is computed by forming the adjoint matrices and taking their product trace.

## The Math

### Lie Groups and Algebras

A **Lie group** G is a smooth manifold with a group structure where multiplication and inversion are smooth. The **Lie algebra** g = T_eG is the tangent space at the identity. For matrix groups, g is a subspace of n×n matrices closed under the commutator bracket.

### The Exponential Map

exp: g → G is defined by exp(X) = Σ X^k/k!. It maps the Lie algebra to the Lie group: for so(n), exp produces rotation matrices. The inverse (where it exists) is the matrix logarithm.

### Baker-Campbell-Hausdorff

BCH gives the product of exponentials as a single exponential: exp(X)exp(Y) = exp(Z) where Z = X + Y + ½[X,Y] + 1/12[X,[X,Y]] − 1/12[Y,[X,Y]] + … This is an infinite series of nested commutators.

### Adjoint Representation

- **Ad(g): G → GL(g)** sends g to the conjugation map X ↦ gXg⁻¹.
- **ad(X): g → End(g)** is the infinitesimal version: ad(X)Y = [X,Y].
- Key relation: **Ad(exp X) = exp(ad X)**.

### Structure Constants

Given a basis {e₁, …, e_d} for g, define c^k_{ij} by [eᵢ, eⱼ] = Σ_k c^k_{ij} eₖ. These completely determine the Lie algebra. For so(3), they're the Levi-Civita symbol ε_{ijk}.

### Killing Form

B(X,Y) = Tr(ad X ∘ ad Y) is the natural symmetric bilinear form on g. It is:
- **Symmetric**: B(X,Y) = B(Y,X)
- **Invariant**: B([Z,X], Y) + B(X, [Z,Y]) = 0
- **Negative definite** for compact semisimple algebras (so(n), n ≥ 3)

Cartan's criterion: g is semisimple iff B is nondegenerate.

### Root Systems

For a semisimple Lie algebra g, choose a Cartan subalgebra h. The **roots** are the nonzero weights of the adjoint representation of h on g. The classical types are:
- **Aₙ**: su(n+1), dim = n(n+2). Roots: eᵢ − eⱼ. n(n+1)/2 positive roots.
- **Bₙ**: so(2n+1), dim = n(2n+1). Roots: ±eᵢ ± eⱼ, ±eᵢ. n² positive roots.
- **Cₙ**: sp(2n), dim = n(2n+1). Roots: ±eᵢ ± eⱼ, ±2eᵢ. n² positive roots.
- **Dₙ**: so(2n), dim = n(2n−1). Roots: ±eᵢ ± eⱼ. n(n−1) positive roots.

The **Cartan matrix** C_{ij} = 2(αᵢ, αⱼ)/(αᵢ, αᵢ) encodes the angles between simple roots.

### Peter-Weyl Theorem

For a compact Lie group G, L²(G) decomposes as ⊕_π V_π ⊗ V_π*, where π runs over all irreducible unitary representations. The **character** χ_π(g) = Tr(ρ_π(g)) is a class function that uniquely determines π. For SU(2), the Weyl character formula gives χ_n(θ) = sin((n+1)θ)/sin(θ).

## License

MIT
