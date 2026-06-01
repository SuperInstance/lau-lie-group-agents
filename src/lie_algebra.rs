//! Lie algebra: vector space with bracket [X,Y] = XY - YX.

use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};

/// A Lie algebra element (an n×n matrix).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LieAlgebraElement {
    pub matrix: DMatrix<f64>,
}

impl LieAlgebraElement {
    /// Create from a matrix.
    pub fn new(matrix: DMatrix<f64>) -> Self {
        Self { matrix }
    }

    /// Zero element.
    pub fn zero(n: usize) -> Self {
        Self { matrix: DMatrix::zeros(n, n) }
    }

    /// Lie bracket [X, Y] = XY - YX.
    pub fn bracket(x: &LieAlgebraElement, y: &LieAlgebraElement) -> LieAlgebraElement {
        LieAlgebraElement {
            matrix: &x.matrix * &y.matrix - &y.matrix * &x.matrix,
        }
    }

    /// Add two elements.
    pub fn add(&self, other: &LieAlgebraElement) -> LieAlgebraElement {
        LieAlgebraElement {
            matrix: &self.matrix + &other.matrix,
        }
    }

    /// Scale by scalar.
    pub fn scale(&self, s: f64) -> LieAlgebraElement {
        LieAlgebraElement {
            matrix: self.matrix.scale(s),
        }
    }

    /// Dimension of the matrix.
    pub fn dim(&self) -> usize {
        self.matrix.nrows()
    }

    /// Check if this is traceless (for sl(n), su(n)).
    pub fn is_traceless(&self) -> bool {
        self.matrix.trace().abs() < 1e-10
    }

    /// Check if this is antisymmetric (for so(n)).
    pub fn is_antisymmetric(&self) -> bool {
        (&self.matrix + &self.matrix.transpose()).norm() < 1e-8
    }

    /// Frobenius norm.
    pub fn norm(&self) -> f64 {
        self.matrix.norm()
    }
}

/// The Lie algebra so(n): antisymmetric n×n matrices.
pub struct SoAlgebra;

impl SoAlgebra {
    /// Basis element for so(n): E_{ij} - E_{ji} for i < j.
    pub fn basis_element(n: usize, i: usize, j: usize) -> LieAlgebraElement {
        let mut m = DMatrix::zeros(n, n);
        m[(i, j)] = 1.0;
        m[(j, i)] = -1.0;
        LieAlgebraElement::new(m)
    }

    /// Full basis for so(n): dimension n(n-1)/2.
    pub fn basis(n: usize) -> Vec<LieAlgebraElement> {
        let mut basis = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                basis.push(Self::basis_element(n, i, j));
            }
        }
        basis
    }

    /// Dimension of so(n).
    pub fn dim(n: usize) -> usize {
        n * (n - 1) / 2
    }
}

/// The Lie algebra sp(2n): matrices X such that XJ + JX^T = 0.
pub struct SpAlgebra;

impl SpAlgebra {
    /// Check membership in sp(2n).
    pub fn is_member(m: &DMatrix<f64>, n: usize) -> bool {
        let dim = 2 * n;
        let j = {
            let mut j = DMatrix::zeros(dim, dim);
            for i in 0..n {
                j[(i, n + i)] = 1.0;
                j[(n + i, i)] = -1.0;
            }
            j
        };
        let check = m * &j + &j * &m.transpose();
        check.norm() < 1e-8
    }
}

/// Verify Jacobi identity: [X,[Y,Z]] + [Y,[Z,X]] + [Z,[X,Y]] = 0.
pub fn verify_jacobi(
    x: &LieAlgebraElement,
    y: &LieAlgebraElement,
    z: &LieAlgebraElement,
) -> bool {
    let yz = LieAlgebraElement::bracket(y, z);
    let zx = LieAlgebraElement::bracket(z, x);
    let xy = LieAlgebraElement::bracket(x, y);

    let term1 = LieAlgebraElement::bracket(x, &yz);
    let term2 = LieAlgebraElement::bracket(y, &zx);
    let term3 = LieAlgebraElement::bracket(z, &xy);

    let sum = term1.add(&term2).add(&term3);
    sum.norm() < 1e-8
}

/// Verify antisymmetry: [X,Y] = -[Y,X].
pub fn verify_antisymmetry(x: &LieAlgebraElement, y: &LieAlgebraElement) -> bool {
    let xy = LieAlgebraElement::bracket(x, y);
    let yx = LieAlgebraElement::bracket(y, x);
    (xy.matrix + yx.matrix).norm() < 1e-8
}

/// Verify bilinearity: [aX+bY, Z] = a[X,Z] + b[Y,Z].
pub fn verify_bilinearity(
    x: &LieAlgebraElement,
    y: &LieAlgebraElement,
    z: &LieAlgebraElement,
    a: f64,
    b: f64,
) -> bool {
    let lhs = LieAlgebraElement::bracket(&x.scale(a).add(&y.scale(b)), z);
    let rhs = LieAlgebraElement::bracket(x, z).scale(a).add(&LieAlgebraElement::bracket(y, z).scale(b));
    (lhs.matrix - rhs.matrix).norm() < 1e-8
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bracket_antisymmetry() {
        let x = LieAlgebraElement::new(DMatrix::from_vec(2, 2, vec![0.0, 1.0, -1.0, 0.0]));
        let y = LieAlgebraElement::new(DMatrix::from_vec(2, 2, vec![0.0, 0.0, 1.0, 0.0]));
        assert!(verify_antisymmetry(&x, &y));
    }

    #[test]
    fn test_jacobi_identity() {
        let x = LieAlgebraElement::new(DMatrix::from_vec(2, 2, vec![0.0, 1.0, -1.0, 0.0]));
        let y = LieAlgebraElement::new(DMatrix::from_vec(2, 2, vec![0.0, 0.0, 1.0, 0.0]));
        let z = LieAlgebraElement::new(DMatrix::from_vec(2, 2, vec![1.0, 0.0, 0.0, -1.0]));
        assert!(verify_jacobi(&x, &y, &z));
    }

    #[test]
    fn test_bilinearity() {
        let x = LieAlgebraElement::new(DMatrix::from_vec(2, 2, vec![0.0, 1.0, -1.0, 0.0]));
        let y = LieAlgebraElement::new(DMatrix::from_vec(2, 2, vec![0.0, 0.0, 1.0, 0.0]));
        let z = LieAlgebraElement::new(DMatrix::from_vec(2, 2, vec![1.0, 0.0, 0.0, -1.0]));
        assert!(verify_bilinearity(&x, &y, &z, 2.0, 3.0));
    }

    #[test]
    fn test_so_basis_dim() {
        assert_eq!(SoAlgebra::dim(3), 3);
        assert_eq!(SoAlgebra::dim(4), 6);
    }

    #[test]
    fn test_so_basis_count() {
        let basis = SoAlgebra::basis(3);
        assert_eq!(basis.len(), 3);
    }

    #[test]
    fn test_so_basis_antisymmetric() {
        let e = SoAlgebra::basis_element(3, 0, 1);
        assert!(e.is_antisymmetric());
    }

    #[test]
    fn test_so_basis_traceless() {
        let e = SoAlgebra::basis_element(3, 0, 1);
        assert!(e.is_traceless());
    }

    #[test]
    fn test_so_bracket() {
        let e01 = SoAlgebra::basis_element(3, 0, 1);
        let e12 = SoAlgebra::basis_element(3, 1, 2);
        let bracket = LieAlgebraElement::bracket(&e01, &e12);
        // [e01, e12] should be proportional to e02
        assert!(bracket.is_antisymmetric());
    }

    #[test]
    fn test_zero_bracket() {
        let x = LieAlgebraElement::zero(3);
        let y = SoAlgebra::basis_element(3, 0, 1);
        let bracket = LieAlgebraElement::bracket(&x, &y);
        assert!(bracket.norm() < 1e-10);
    }

    #[test]
    fn test_self_bracket_zero() {
        let x = SoAlgebra::basis_element(3, 0, 1);
        let bracket = LieAlgebraElement::bracket(&x, &x);
        assert!(bracket.norm() < 1e-10);
    }

    #[test]
    fn test_scale() {
        let x = SoAlgebra::basis_element(3, 0, 1);
        let scaled = x.scale(2.0);
        assert!((scaled.norm() - 2.0 * x.norm()).abs() < 1e-10);
    }

    #[test]
    fn test_sp_algebra_membership() {
        let n = 1;
        let dim = 2;
        let m = DMatrix::from_vec(dim, dim, vec![0.0, 1.0, -1.0, 0.0]);
        assert!(SpAlgebra::is_member(&m, n));
    }

    #[test]
    fn test_sp_algebra_not_member() {
        let m = DMatrix::identity(2, 2);
        assert!(!SpAlgebra::is_member(&m, 1));
    }

    #[test]
    fn test_jacobi_so3() {
        let e0 = SoAlgebra::basis_element(3, 0, 1);
        let e1 = SoAlgebra::basis_element(3, 1, 2);
        let e2 = SoAlgebra::basis_element(3, 0, 2);
        assert!(verify_jacobi(&e0, &e1, &e2));
    }

    #[test]
    fn test_add() {
        let e0 = SoAlgebra::basis_element(3, 0, 1);
        let e1 = SoAlgebra::basis_element(3, 1, 2);
        let sum = e0.add(&e1);
        assert_eq!(sum.dim(), 3);
    }

    #[test]
    fn test_serde() {
        let e = SoAlgebra::basis_element(3, 0, 1);
        let json = serde_json::to_string(&e).unwrap();
        let back: LieAlgebraElement = serde_json::from_str(&json).unwrap();
        assert_eq!(back.dim(), 3);
    }
}
