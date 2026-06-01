# lau-lie-group-agents

**Lie group theory for agents: SO(n), SU(n), GL(n), Sp(2n), SE(3), Lie algebras, exponential map, BCH formula, adjoint representation, structure constants, Killing form, root systems, and Peter-Weyl theorem.**

A Rust library implementing the core structures of Lie theory — from classical matrix groups and their algebras to representation theory. Lie groups provide the mathematical language for continuous symmetry, making them foundational for physics, robotics, control theory, and geometric deep learning.

## What This Does

This library provides 9 modules covering Lie theory from foundations to representation theory:

1. **Lie groups** — Classical matrix groups: SO(n), SU(n), GL(n), Sp(2n), SE(3) with multiplication, inversion, and membership tests.
2. **Lie algebras** — The tangent space at identity with Lie bracket [X,Y] = XY−YX, basis construction for so(n) and sp(2n), and axiomatic verification (Jacobi, antisymmetry, bilinearity).
3. **Exponential map** — The bridge exp: 𝔤 → G from algebra to group via matrix exponential (scaling-and-squaring) and logarithm (inverse scaling-and-squaring).
4. **Baker-Campbell-Hausdorff (BCH)** — The formula log(exp(X)exp(Y)) = X + Y + ½[X,Y] + 1/12[X,[X,Y]] − 1/12[Y,[X,Y]] + ... at orders 1–3.
5. **Adjoint representation** — Ad(g)X = gXg⁻¹ (group action on algebra) and ad(X)Y = [X,Y] (infinitesimal), with matrix form and composition law Ad(gh) = Ad(g)Ad(h).
6. **Structure constants** — The tensor cᵏᵢⱼ encoding [eᵢ,eⱼ] = Σ cᵏᵢⱼ eₖ, computed automatically from any basis, with Jacobi and antisymmetry verification.
7. **Killing form** — The invariant bilinear form B(X,Y) = Tr(ad(X)ad(Y)), with negative-definiteness check (Cartan's criterion for semisimplicity).
8. **Root systems** — Classical types Aₙ (su(n+1)), Bₙ (so(2n+1)), Cₙ (sp(2n)), Dₙ (so(2n)) with simple roots, positive roots, and Cartan matrices.
9. **Peter-Weyl theorem** — Decomposition of L²(G) into irreducible representations, characters, Weyl character formula for SU(2), and Weyl dimension formula.

## Key Idea

A **Lie group** G is a group that is also a smooth manifold — think rotation matrices (SO(n)) or invertible matrices (GL(n)). Its **Lie algebra** 𝔤 is the tangent space at the identity, equipped with the Lie bracket [X,Y] = XY − YX.

The magic: the entire local structure of G is encoded in 𝔤. The exponential map bridges them — exp(X) sends algebra elements to group elements, and the BCH formula tells you how to multiply in the algebra: exp(X)exp(Y) = exp(BCH(X,Y)).

This library makes that bridge explicit and computational. You can:
- Construct group elements and multiply them
- Work in the algebra with brackets and structure constants
- Map between group and algebra via exp/log
- Verify algebraic identities (Jacobi, antisymmetry, invariance)
- Classify algebras via root systems and the Killing form

## Install

```toml
[dependencies]
lau-lie-group-agents = "0.1.0"
```

Or:

```sh
cargo add lau-lie-group-agents
```

Requires Rust 2021 edition. Dependencies: `nalgebra` (with serde), `serde`/`serde_json`, `num-complex`.

## Quick Start

### Classical Lie Groups

```rust
use lau_lie_group_agents::lie_groups::{SO, GL, SU, Sp, SE3, LieGroupElement, GroupType};
use nalgebra::DVector;

// SO(3) rotations
let rx = SO::rotx(std::f64::consts::PI / 2.0); // 90° around x
let ry = SO::roty(1.0);                         // 1 rad around y
let rz = SO::rotz(0.5);                         // 0.5 rad around z
assert!(SO::is_member(&rx.matrix));              // Verify membership

// Compose: SO(2) rotations add angles
let r1 = SO::so2(1.0);
let r2 = SO::so2(2.0);
let r3 = r1.multiply(&r2);
let r_direct = SO::so2(3.0);
assert!((r3.matrix - r_direct.matrix).norm() < 1e-10);

// Inverse
let inv = rx.inverse().unwrap();
let product = rx.multiply(&inv); // ≈ identity

// SE(3) rigid body transforms
let t = SE3::translation(1.0, 2.0, 3.0);
let p = DVector::from_vec(vec![0.0, 0.0, 0.0]);
let moved = SE3::transform_point(&t.matrix, &p);
// moved = (1, 2, 3)

// Sp(2n) symplectic matrices
let sp = Sp::symplectic_rotation(2, 0.5);
assert!(Sp::is_member(&sp.matrix, 2));
```

### Lie Algebras

```rust
use lau_lie_group_agents::lie_algebra::{LieAlgebraElement, SoAlgebra, SpAlgebra,
    verify_jacobi, verify_antisymmetry, verify_bilinearity};

// Basis for so(3): 3 antisymmetric matrices (dimension 3)
let basis = SoAlgebra::basis(3);
assert_eq!(basis.len(), 3);

// Lie bracket [X,Y] = XY - YX
let e01 = SoAlgebra::basis_element(3, 0, 1);
let e12 = SoAlgebra::basis_element(3, 1, 2);
let bracket = LieAlgebraElement::bracket(&e01, &e12);
assert!(bracket.is_antisymmetric());

// Verify algebra axioms
let e02 = SoAlgebra::basis_element(3, 0, 2);
assert!(verify_jacobi(&e01, &e12, &e02));       // [X,[Y,Z]] + cyclic = 0
assert!(verify_antisymmetry(&e01, &e12));         // [X,Y] = -[Y,X]
assert!(verify_bilinearity(&e01, &e12, &e02, 2.0, 3.0));

// Sp algebra membership
use nalgebra::DMatrix;
let m = DMatrix::from_vec(2, 2, vec![0.0, 1.0, -1.0, 0.0]);
assert!(SpAlgebra::is_member(&m, 1));
```

### Exponential Map

```rust
use lau_lie_group_agents::exponential_map::{ExponentialMap, matrix_exp, matrix_log};
use lau_lie_group_agents::lie_algebra::SoAlgebra;

// exp maps algebra → group
let map = ExponentialMap::new(3);
let x = SoAlgebra::basis_element(3, 0, 1);
let g = map.apply(&x);
assert!(crate::lie_groups::SO::is_member(&g.matrix)); // exp(so(n)) ⊂ SO(n)

// exp(0) = I
assert!(map.verify_identity());

// Full SO(n) mapping verified
assert!(map.verify_so_mapping());

// Full rotation: exp(2π·X) = I for so(2)
let x2 = SoAlgebra::basis_element(2, 0, 1);
let full_rot = matrix_exp(&x2.matrix.scale(2.0 * std::f64::consts::PI));
assert!((full_rot - DMatrix::identity(2, 2)).norm() < 1e-8);

// Log roundtrip: exp(log(M)) = M
let m = DMatrix::from_vec(2, 2, vec![2.0, 0.0, 0.0, 3.0]);
let exp_log = matrix_exp(&matrix_log(&m));
assert!((exp_log - m).norm() < 1e-6);
```

### Baker-Campbell-Hausdorff

```rust
use lau_lie_group_agents::bch::{bch, verify_bch, BCHOrder};
use lau_lie_group_agents::lie_algebra::LieAlgebraElement;
use nalgebra::DMatrix;

let x = LieAlgebraElement::new(DMatrix::from_vec(2, 2, vec![0.0, 0.1, -0.1, 0.0]));
let y = LieAlgebraElement::new(DMatrix::from_vec(2, 2, vec![0.0, 0.0, 0.1, 0.0]));

// BCH at different orders (higher = more accurate)
let first  = bch(&x, &y, BCHOrder::First);   // X + Y
let second = bch(&x, &y, BCHOrder::Second);  // X + Y + ½[X,Y]
let third  = bch(&x, &y, BCHOrder::Third);   // + 1/12[X,[X,Y]] - 1/12[Y,[X,Y]]

// Verify: exp(X)exp(Y) ≈ exp(BCH(X,Y))
assert!(verify_bch(&x, &y, BCHOrder::Third));

// Commutative case: [X,Y]=0 → BCH = X+Y exactly
let zero = LieAlgebraElement::zero(2);
let result = bch(&zero, &y, BCHOrder::Third);
assert!((result.matrix - y.matrix).norm() < 1e-10);
```

### Adjoint Representation

```rust
use lau_lie_group_agents::adjoint::AdjointRep;
use lau_lie_group_agents::lie_groups::SO;
use lau_lie_group_agents::lie_algebra::SoAlgebra;

// Ad(g)X = gXg⁻¹: conjugation action
let g = SO::rotz(0.5);
let x = SoAlgebra::basis_element(3, 0, 1);
let ad_g_x = AdjointRep::ad_g(&g, &x);
assert!(ad_g_x.is_antisymmetric()); // stays in algebra

// ad(X)Y = [X,Y]: infinitesimal adjoint
let y = SoAlgebra::basis_element(3, 1, 2);
let ad_x_y = AdjointRep::ad_x(&x, &y);

// Matrix representation of ad(X) in so(3) basis
let basis = SoAlgebra::basis(3);
let ad_matrix = AdjointRep::ad_matrix(&x, &basis);
assert_eq!(ad_matrix.shape(), (3, 3));

// Ad(gh) = Ad(g)∘Ad(h) (composition)
let g1 = SO::rotz(0.3);
let g2 = SO::roty(0.5);
let g12 = g1.multiply(&g2);
let ad_g12 = AdjointRep::ad_g(&g12, &x);
let ad_seq = AdjointRep::ad_g(&g1, &AdjointRep::ad_g(&g2, &x));
assert!((ad_g12.matrix - ad_seq.matrix).norm() < 1e-8);

// Trace of ad = 0 for semisimple algebras
assert!(AdjointRep::trace_ad(&x, &basis).abs() < 1e-8);
```

### Structure Constants

```rust
use lau_lie_group_agents::structure_constants::StructureConstants;

// Structure constants for so(3): [eᵢ,eⱼ] = Σ cᵏᵢⱼ eₖ
let sc = StructureConstants::for_so(3);
assert_eq!(sc.dim(), 3);

// Verify algebraic properties
assert!(sc.verify_antisymmetry()); // cᵏᵢⱼ = -cᵏⱼᵢ
assert!(sc.verify_jacobi());       // Jacobi from constants

// For so(3), cᵏᵢⱼ = εᵢⱼₖ (Levi-Civita)
// [e₁,e₂] = e₃ (up to sign depending on normalization)
let c = sc.get(0, 1, 2);
assert!((c.abs() - 1.0).abs() < 1e-8);

// so(2) is abelian: all constants are 0
let sc2 = StructureConstants::for_so(2);
assert!(sc2.constants[0][0][0].abs() < 1e-10);

// so(4) has dimension 6
let sc4 = StructureConstants::for_so(4);
assert_eq!(sc4.dim(), 6);
assert!(sc4.verify_jacobi());
```

### Killing Form

```rust
use lau_lie_group_agents::killing_form::KillingForm;
use lau_lie_group_agents::lie_algebra::SoAlgebra;

let kf = KillingForm::for_so(3);
let basis = SoAlgebra::basis(3);

// B(X,Y) = Tr(ad(X)ad(Y))
let b_matrix = kf.matrix(&basis);
assert_eq!(b_matrix.shape(), (3, 3));

// so(n) is semisimple: Killing form is negative definite
assert!(kf.is_negative_definite(3));
assert!(kf.is_negative_definite(4));

// Invariant properties
assert!(kf.verify_bilinearity(&basis));
assert!(kf.verify_symmetry(&basis));      // B(X,Y) = B(Y,X)
assert!(kf.verify_invariance(&basis));     // B([Z,X],Y) + B(X,[Z,Y]) = 0

// Diagonal is negative (for compact semisimple groups)
for i in 0..3 {
    assert!(b_matrix[(i, i)] < 0.0);
}
```

### Root Systems

```rust
use lau_lie_group_agents::root_systems::{RootSystem, RootSystemType, count_positive_roots};

// A₂ (su(3)): 2 simple roots, 3 positive roots, 6 total
let a2 = RootSystem::type_a(2);
assert_eq!(a2.num_simple_roots(), 2);
assert_eq!(a2.num_positive_roots(), 3);
assert_eq!(a2.num_roots(), 6);

// B₃ (so(7)): 9 positive roots
let b3 = RootSystem::type_b(3);
assert_eq!(b3.num_positive_roots(), 9);

// D₄ (so(8)): 12 positive roots
let d4 = RootSystem::type_d(4);
assert_eq!(d4.num_positive_roots(), 12);

// Cartan matrix: 2's on diagonal
let c = a2.cartan_matrix();
assert!((c[(0, 0)] - 2.0).abs() < 1e-8);
assert!(a2.verify_cartan_diagonal());

// Count positive roots by formula
assert_eq!(count_positive_roots(RootSystemType::A, 3), 6);  // n(n+1)/2
assert_eq!(count_positive_roots(RootSystemType::B, 3), 9);  // n²
assert_eq!(count_positive_roots(RootSystemType::D, 4), 12); // n(n-1)
```

### Peter-Weyl Theorem

```rust
use lau_lie_group_agents::peter_weyl::{
    Representation, PeterWeylDecomposition, su2_character, su_n_dimension
};

// Decompose L²(SO(3)) into irreps
let pw = PeterWeylDecomposition::for_so(3);
assert_eq!(pw.num_representations(), 3);     // trivial, standard, adjoint
assert_eq!(pw.total_dimension(), 19);         // 1² + 3² + 3²

// Character of a representation: χ(g) = Tr(ρ(g))
let rep = Representation::standard(2);
let g = DMatrix::from_vec(2, 2, vec![
    theta.cos(), theta.sin(), -theta.sin(), theta.cos(),
]);
let chi = rep.character(&g); // = 2cos(θ) for SO(2)

// Weyl character formula for SU(2): χₙ(θ) = sin((n+1)θ)/sin(θ)
assert!((su2_character(0, std::f64::consts::PI/4.0) - 1.0).abs() < 1e-10);
assert!((su2_character(1, 0.0) - 2.0).abs() < 1e-10); // dim of spin-½ rep

// Matrix coefficient orthogonality
let rep1 = Representation::standard(2);
let rep2 = Representation::standard(2);
let inner = PeterWeylDecomposition::matrix_coefficient_inner_product(
    &g, &rep1, &rep2, 0, 0, 0, 0,
);
assert!((inner - 0.5).abs() < 1e-10); // 1/dim
```

## API Reference

### `lie_groups`

| Type | Description |
|------|-------------|
| `LieGroupElement` | Matrix element of a Lie group |
| `GroupType` | Enum: SO(n), SU(n), GL(n), Sp(n), SE3 |
| `SO` | Special Orthogonal: rotation matrices |
| `GL` | General Linear: invertible matrices |
| `SU` | Special Unitary (realified) |
| `Sp` | Symplectic: preserves ω = J |
| `SE3` | Special Euclidean: rigid transforms |

### `lie_algebra`

| Type | Description |
|------|-------------|
| `LieAlgebraElement` | Matrix in a Lie algebra |
| `SoAlgebra` | so(n): antisymmetric matrices |
| `SpAlgebra` | sp(2n): XJ + JX^T = 0 |
| `verify_jacobi` | Jacobi identity check |
| `verify_antisymmetry` | [X,Y] = −[Y,X] check |
| `verify_bilinearity` | Bracket bilinearity check |

### `exponential_map`

| Function | Description |
|----------|-------------|
| `matrix_exp` | Matrix exponential (scaling + squaring) |
| `matrix_log` | Matrix logarithm (inverse scaling + squaring) |
| `ExponentialMap` | Wrapper: apply, verify_identity, verify_so_mapping, exp_t |

### `bch`

| Function | Description |
|----------|-------------|
| `bch(X, Y, order)` | BCH formula at order 1/2/3 |
| `verify_bch(X, Y, order)` | Verify exp(X)exp(Y) ≈ exp(BCH) |
| `bch_series(X, Y, depth)` | Commutator series to arbitrary depth |
| `BCHOrder` | First, Second, Third |

### `adjoint`

| Function | Description |
|----------|-------------|
| `AdjointRep::ad_g(g, X)` | Ad(g)X = gXg⁻¹ |
| `AdjointRep::ad_x(X, Y)` | ad(X)Y = [X,Y] |
| `AdjointRep::ad_matrix(X, basis)` | Matrix of ad(X) in basis |
| `AdjointRep::trace_ad(X, basis)` | Tr(ad(X)), zero for semisimple |

### `structure_constants`

| Type | Description |
|------|-------------|
| `StructureConstants` | Tensor cᵏᵢⱼ for [eᵢ,eⱼ] = Σ cᵏᵢⱼ eₖ |

Key methods: `for_so(n)`, `from_basis`, `get(i,j,k)`, `verify_antisymmetry`, `verify_jacobi`.

### `killing_form`

| Type | Description |
|------|-------------|
| `KillingForm` | B(X,Y) = Tr(ad(X)ad(Y)) |

Key methods: `for_so(n)`, `evaluate(X,Y,basis)`, `matrix(basis)`, `is_negative_definite`, `verify_bilinearity`, `verify_symmetry`, `verify_invariance`.

### `root_systems`

| Type | Description |
|------|-------------|
| `RootSystem` | Classical root system with positive/simple roots |
| `RootSystemType` | A, B, C, D |
| `count_positive_roots` | Formula: n(n+1)/2, n², n², n(n−1) |

Key methods: `type_a/b/c_d(n)`, `cartan_matrix`, `verify_cartan_diagonal`, `num_roots`.

### `peter_weyl`

| Type | Description |
|------|-------------|
| `Representation` | Unitary rep with dimension and type |
| `RepType` | Trivial, Standard, Adjoint, Fundamental, Custom |
| `PeterWeylDecomposition` | L²(G) = ⊕ Vₐ ⊗ Vₐ* |
| `su2_character(n, θ)` | Weyl character formula for SU(2) |
| `su_n_dimension(n, weights)` | Weyl dimension formula |

## How It Works

### Lie Groups as Matrix Groups

Each classical group is represented as a matrix type with group operations (multiplication = matrix multiply, inverse = matrix inverse). Membership tests verify the defining property:
- **SO(n)**: M^TM = I, det(M) = 1
- **Sp(2n)**: M^TJM = J (symplectic condition)
- **GL(n)**: det(M) ≠ 0

### Lie Bracket and Axioms

The bracket [X,Y] = XY − YX satisfies three axioms verified numerically:
1. **Antisymmetry**: [X,Y] = −[Y,X]
2. **Bilinearity**: [aX+bY, Z] = a[X,Z] + b[Y,Z]
3. **Jacobi identity**: [X,[Y,Z]] + [Y,[Z,X]] + [Z,[X,Y]] = 0

### Exponential Map

The matrix exponential exp(X) = I + X + X²/2! + X³/3! + ... is computed via scaling-and-squaring: scale X until ‖X‖ < 1, compute the Taylor series, then square the result back. The logarithm uses inverse scaling-and-squaring.

Key property: exp maps so(n) → SO(n). Any antisymmetric matrix exponentiates to an orthogonal matrix.

### BCH Formula

For small X, Y near 0: log(exp(X)exp(Y)) can be expressed entirely in terms of X, Y, and nested commutators. The library computes this to third order:
- Z₁ = X + Y
- Z₂ = X + Y + ½[X,Y]
- Z₃ = X + Y + ½[X,Y] + 1/12[X,[X,Y]] − 1/12[Y,[X,Y]]

Higher orders converge for ‖X‖, ‖Y‖ < π.

### Adjoint Representation

The adjoint representation Ad: G → GL(𝔤) maps group elements to linear maps on the algebra via conjugation: Ad(g)X = gXg⁻¹. Its derivative ad(X)Y = [X,Y] is the bracket itself. The key identity: Ad(exp(X)) = exp(ad(X)).

### Structure Constants

For a basis {eᵢ} of 𝔤, the bracket decomposes as [eᵢ,eⱼ] = Σₖ cᵏᵢⱼ eₖ. The constants cᵏᵢⱼ encode the entire algebra. For so(3), they are the Levi-Civita symbol εᵢⱼₖ.

### Killing Form

B(X,Y) = Tr(ad(X)ad(Y)) is an invariant symmetric bilinear form on the algebra. **Cartan's criterion**: 𝔤 is semisimple iff B is nondegenerate. For compact semisimple algebras (like so(n)), B is negative definite.

### Root Systems

For a semisimple algebra with Cartan subalgebra 𝔥, roots are the nonzero weights of the adjoint action of 𝔥. Classical types:
- **Aₙ**: su(n+1), positive roots = {eᵢ − eⱼ : i < j}, count = n(n+1)/2
- **Bₙ**: so(2n+1), adds short roots {eᵢ}, count = n²
- **Cₙ**: sp(2n), adds long roots {2eᵢ}, count = n²
- **Dₙ**: so(2n), count = n(n−1)

The Cartan matrix Aᵢⱼ = 2⟨αᵢ,αⱼ⟩/⟨αᵢ,αᵢ⟩ has 2's on diagonal and non-positive integers off-diagonal.

### Peter-Weyl Theorem

For a compact Lie group G, L²(G) decomposes as ⊕_π Vₐ ⊗ Vₐ* where π runs over irreducible unitary representations. Matrix coefficients of different irreps are orthogonal in L²(G). The Weyl character formula gives χₙ(θ) = sin((n+1)θ)/sin(θ) for SU(2).

## The Math

### Lie Group → Lie Algebra

At the identity e ∈ G, the tangent space TₑG = 𝔤 is the Lie algebra. For matrix groups, algebra elements are matrices: so(n) = {X : X^T = −X}, sp(2n) = {X : XJ + JX^T = 0}.

### Exponential Map

$$\exp(X) = \sum_{k=0}^{\infty} \frac{X^k}{k!}$$

Inverse: $\log(M) = (M-I) - \frac{(M-I)^2}{2} + \frac{(M-I)^3}{3} - \cdots$ near identity.

### BCH Formula

$$\log(\exp(X)\exp(Y)) = X + Y + \frac{1}{2}[X,Y] + \frac{1}{12}[X,[X,Y]] - \frac{1}{12}[Y,[X,Y]] + \cdots$$

### Structure Constants

$$[e_i, e_j] = \sum_k c_{ij}^k \, e_k$$

Jacobi in terms of constants: $\sum_m (c_{ij}^m c_{mk}^l + c_{jk}^m c_{mi}^l + c_{ki}^m c_{mj}^l) = 0$.

### Killing Form

$$B(X, Y) = \mathrm{Tr}(\mathrm{ad}(X) \circ \mathrm{ad}(Y))$$

Properties: symmetric, bilinear, invariant ($B([Z,X],Y) + B(X,[Z,Y]) = 0$). Negative definite ⟺ compact semisimple.

### Root Systems

For type Aₙ: simple roots αᵢ = eᵢ − eᵢ₊₁, Cartan matrix has 2 on diagonal, −1 on super/sub-diagonal.

Positive roots count:
| Type | Group | Positive roots |
|------|-------|---------------|
| Aₙ | su(n+1) | n(n+1)/2 |
| Bₙ | so(2n+1) | n² |
| Cₙ | sp(2n) | n² |
| Dₙ | so(2n) | n(n−1) |

### Peter-Weyl

$$L^2(G) \cong \bigoplus_{\pi \in \hat{G}} V_\pi \otimes V_\pi^*$$

Orthogonality: $\langle \rho_1(g)_{ij}, \rho_2(g)_{kl} \rangle = \frac{\delta_{\pi_1\pi_2}\delta_{ik}\delta_{jl}}{\dim \pi}$

SU(2) Weyl character: $\chi_n(\theta) = \frac{\sin((n+1)\theta)}{\sin(\theta)}$

## Test Coverage

118 tests across 9 modules:
- **lie_groups** (18): SO identity/rotations/composition/inverse, GL identity/diagonal/membership, SU(2), Sp identity/J/rotation, SE3 identity/translation/transform/compose, GroupType dimensions, Serde
- **lie_algebra** (16): bracket antisymmetry/Jacobi/bilinearity, so(n) basis dim/count/antisymmetry/traceless/bracket, zero/self bracket, scale, sp membership, add, Serde
- **exponential_map** (10): exp(0), exp(I), exp(antisymmetric)→SO, exp so(3), identity verification, SO mapping, exp_t, full 2π rotation, log roundtrip, Serde
- **bch** (10): orders 1/2/3, commutative case, verification 1/2/3, series, higher order accuracy, Serde
- **adjoint** (10): Ad(g)/ad(X) basic, identity, matrix, trace=0, preserves algebra, Jacobi identity, composition, so(2) trivial, Serde
- **structure_constants** (11): so(3) construction/antisymmetry/Jacobi, so(2) abelian, so(4) construction/Jacobi/antisymmetry, Levi-Civita, nonzero bracket, custom basis, Serde
- **killing_form** (12): so(3) construction/negative-definite/bilinearity/symmetry/invariance, so(2), so(4) construction/negative-definite/symmetry, matrix, diagonal negativity, Serde
- **root_systems** (15): A₁/A₂/A₃/B₂/C₃/D₃/D₄ roots, Cartan matrix, diagonal verification, counts A/B/C/D, total roots, Serde
- **peter_weyl** (16): trivial/standard/adjoint reps, character identity/rotation, SO(3)/SO(2) decomposition, orthogonality, total dimension, SU(2) character (spin 0/½/1), matrix coefficient orth/same, SU(n) dimension, Serde

## License

MIT
