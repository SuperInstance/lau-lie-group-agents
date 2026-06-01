# lau-lie-group-agents

> Lie groups, Lie algebras, exponential maps, BCH formula, root systems, and the Peter-Weyl theorem for agent symmetry

## What This Does

This crate implements Lie group theory — the mathematics of continuous symmetry — for agent systems. It covers classical Lie groups (SO(n), SU(n), GL(n), Sp(2n), SE(3)), Lie algebras with bracket operations, the exponential map (algebra → group), Baker-Campbell-Hausdorff formula, adjoint representation, structure constants, Killing form, root systems (A_n, B_n, C_n, D_n), and the Peter-Weyl theorem on irreducible representations.

## The Key Idea

Symmetry is the deepest idea in physics and mathematics. A Lie group is a group that's also a smooth manifold — rotations, Lorentz transforms, gauge transformations are all Lie groups. The Lie algebra (the tangent space at the identity) captures the group locally via the exponential map. BCH tells you how to multiply in the algebra: log(exp(X)exp(Y)) = X + Y + ½[X,Y] + ... For agents, this means: if an agent's state transforms under a symmetry group, the Lie algebra tells you how small perturbations evolve.

## Install

```toml
[dependencies]
lau-lie-group-agents = { git = "https://github.com/SuperInstance/lau-lie-group-agents" }
```

## Quick Start

```rust
use lau_lie_group_agents::*;
use nalgebra::DMatrix;

// === Lie Group: SO(3) rotations ===
let g1 = LieGroupElement::new(
    DMatrix::from_row_slice(3, 3, &[1,0,0, 0,1,0, 0,0,1]),
    GroupType::SO(3),
);
let g2 = g1.multiply(&g1);
let inv = g1.inverse().unwrap();
assert_eq!(g1.multiply(&inv), LieGroupElement::identity(GroupType::SO(3)));

// === Lie Algebra: so(3) ===
let X = LieAlgebraElement::new(DMatrix::from_row_slice(3, 3,
    &[0,-1,0, 1,0,0, 0,0,0])); // Generator of rotation around z
let Y = LieAlgebraElement::new(DMatrix::from_row_slice(3, 3,
    &[0,0,1, 0,0,0, -1,0,0])); // Generator of rotation around y

// Lie bracket [X, Y] = XY - YX
let bracket = LieAlgebraElement::bracket(&X, &Y);
println!("[X,Y] = {:?}", bracket.matrix);

// === Exponential Map ===
let exp_X = exponential_map::matrix_exp(&X.matrix);
println!("exp(X) = rotation matrix:\n{:?}", exp_X);

// === BCH Formula ===
let bch_result = bch::bch(&X, &Y, bch::BCHOrder::Third);
println!("BCH(X,Y) ≈ {:?}", bch_result.matrix);

// === Adjoint Representation ===
let basis = lie_algebra::SoAlgebra::basis(3);
let ad_mat = AdjointRep::ad_matrix(&X, &basis);
println!("ad(X) matrix:\n{:?}", ad_mat);

// === Structure Constants ===
let sc = StructureConstants::from_basis(&basis);
println!("c^k_{{12}} = {}", sc.constants[0][1][2]);

// === Killing Form ===
let kf = KillingForm::for_so(3);
let bilinear = kf.evaluate(&X, &Y, &basis);
println!("B(X,Y) = {}", bilinear);
println!("Negative definite (semisimple): {}", kf.is_negative_definite(3));

// === Root System A_2 (su(3)) ===
let roots = RootSystem::type_a(2);
println!("Simple roots: {:?}", roots.simple_roots);
println!("Positive roots: {} total", roots.positive_roots.len());

// === Peter-Weyl: representations ===
let trivial = Representation::trivial();
let standard = Representation::standard(3);
let adjoint = Representation::adjoint(3);
println!("χ(g) = {}", standard.character(&exp_X));
println!("Irreducible (heuristic): {}", standard.is_irreducible_heuristic());
```

## API Reference

### `lie_groups`

| Type | Description |
|------|-------------|
| `LieGroupElement::new(matrix, group_type)` | A group element. |
| `LieGroupElement::identity(group_type)` | Identity element. |
| `multiply(other)` | Group multiplication. |
| `inverse()` | Group inverse. |
| `GroupType` | Enum: `SO(n)`, `SU(n)`, `GL(n)`, `Sp(n)`, `SE3`. |

### `lie_algebra`

| Type | Description |
|------|-------------|
| `LieAlgebraElement::new(matrix)` | An algebra element (n×n matrix). |
| `bracket(x, y)` | [X,Y] = XY − YX. |
| `add(other)`, `scale(s)` | Vector space operations. |
| `is_antisymmetric()` | For so(n). |
| `is_traceless()` | For sl(n), su(n). |
| `norm()` | Frobenius norm. |
| `SoAlgebra::basis(n)` | Standard basis for so(n). |
| `SoAlgebra::dim(n)` | Dimension n(n−1)/2. |

### `exponential_map`

| Function | Description |
|----------|-------------|
| `matrix_exp(m)` | Compute exp(M) via scaling-and-squaring. |
| `matrix_log(m)` | Compute log(M) (inverse). |

### `bch`

| Function | Description |
|----------|-------------|
| `bch(x, y, order)` | BCH: log(exp(X)exp(Y)) to first/second/third order. |
| `verify_bch(x, y, order)` | Verify exp(X)exp(Y) ≈ exp(BCH(X,Y)). |

### `adjoint`

| Type | Description |
|------|-------------|
| `AdjointRep::ad_g(g, x)` | Ad(g)X = gXg⁻¹. |
| `AdjointRep::ad_x(x, y)` | ad(X)Y = [X,Y]. |
| `AdjointRep::ad_matrix(x, basis)` | Matrix representation of ad(X). |
| `verify_ad_exp(x)` | Check Ad(exp(X)) = exp(ad(X)). |

### `structure_constants`

| Type | Description |
|------|-------------|
| `StructureConstants::from_basis(basis)` | Compute c^k_{ij} where [eᵢ,eⱼ] = Σ c^k_{ij} eₖ. |

### `killing_form`

| Type | Description |
|------|-------------|
| `KillingForm::for_so(n)` | B(X,Y) = Tr(ad(X)ad(Y)). |
| `evaluate(x, y, basis)` | Compute B(X,Y). |
| `matrix(basis)` | Full Killing form matrix. |
| `is_negative_definite(n)` | Check semisimplicity. |

### `root_systems`

| Type | Description |
|------|-------------|
| `RootSystem::type_a(n)` | Aₙ roots (su(n+1)). |
| `RootSystem::type_b(n)` | Bₙ roots (so(2n+1)). |
| `RootSystem::type_c(n)` | Cₙ roots (sp(2n)). |
| `RootSystem::type_d(n)` | Dₙ roots (so(2n)). |
| `simple_roots`, `positive_roots` | Root data. |
| `cartan_matrix()` | Cartan matrix from simple roots. |
| `weyl_group_size()` | Order of the Weyl group. |

### `peter_weyl`

| Type | Description |
|------|-------------|
| `Representation` | Unitary rep with `dim`, `rep_type`. |
| `RepType` | Trivial, Standard, Adjoint, Fundamental, Custom. |
| `character(g)` | χ(g) = Tr(ρ(g)). |
| `is_irreducible_heuristic()` | Check irreducibility. |

## How It Works

1. **Group Elements**: Stored as matrices. Multiplication = matrix multiply. Inverse = matrix inverse.
2. **Algebra Elements**: Also matrices, but with constraints (antisymmetric for so(n), traceless+antihermitian for su(n)).
3. **Exponential Map**: Scaling-and-squaring with Taylor series. Maps algebra elements to group elements.
4. **BCH**: Approximate formula for combining algebra elements: log(e^X e^Y) = X + Y + ½[X,Y] + ...
5. **Adjoint**: The derivative of conjugation. Ad(g)X = gXg⁻¹ at the group level, ad(X)Y = [X,Y] at the algebra level.
6. **Killing Form**: The natural inner product on a Lie algebra. Negative definite ↔ semisimple algebra.
7. **Root Systems**: Encode the internal structure of semisimple Lie algebras. Simple roots generate all roots via reflections.

## The Math

- **Lie Group**: Smooth manifold G with group operation (g,h) ↦ gh that's smooth.
- **Lie Algebra**: Tangent space g = T_eG with bracket [X,Y] = XY − YX satisfying Jacobi identity.
- **Exponential Map**: exp: g → G, exp(X) = Σ Xⁿ/n!.
- **BCH**: log(e^X e^Y) = X + Y + ½[X,Y] + 1/12[X,[X,Y]] − 1/12[Y,[X,Y]] + ...
- **Killing Form**: B(X,Y) = Tr(ad(X) ∘ ad(Y)). Non-degenerate ↔ semisimple.
- **Peter-Weyl**: L²(G) = ⊕_π dim(V_π) · V_π for compact G (generalized Fourier analysis).

## Testing

118 tests covering:
- Group element multiplication and inversion
- Lie bracket Jacobi identity
- Exponential map accuracy (exp(X) ∈ G)
- BCH convergence at each order
- Adjoint representation consistency
- Structure constants antisymmetry
- Killing form bilinearity and negative definiteness for so(n)
- Root system properties (number of roots, Cartan matrix, Weyl group)
- Peter-Weyl representation character computation

## License

MIT
