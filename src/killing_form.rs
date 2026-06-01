//! Killing form: B(X,Y) = Tr(ad(X) ∘ ad(Y)).

use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};

use crate::adjoint::AdjointRep;
use crate::lie_algebra::{LieAlgebraElement, SoAlgebra};

/// The Killing form B(X,Y) = Tr(ad(X)ad(Y)).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KillingForm {
    pub algebra_dim: usize,
    pub matrix_dim: usize,
}

impl KillingForm {
    /// Create for so(n).
    pub fn for_so(n: usize) -> Self {
        Self {
            algebra_dim: SoAlgebra::dim(n),
            matrix_dim: n,
        }
    }

    /// Evaluate B(X,Y) = Tr(ad(X) ∘ ad(Y)).
    pub fn evaluate(&self, x: &LieAlgebraElement, y: &LieAlgebraElement, basis: &[LieAlgebraElement]) -> f64 {
        let ad_x = AdjointRep::ad_matrix(x, basis);
        let ad_y = AdjointRep::ad_matrix(y, basis);
        let product = &ad_x * &ad_y;
        product.trace()
    }

    /// Compute the Killing form matrix B_{ij} = B(eᵢ, eⱼ).
    pub fn matrix(&self, basis: &[LieAlgebraElement]) -> DMatrix<f64> {
        let dim = basis.len();
        let mut b = DMatrix::zeros(dim, dim);
        for i in 0..dim {
            for j in 0..dim {
                b[(i, j)] = self.evaluate(&basis[i], &basis[j], basis);
            }
        }
        b
    }

    /// Check if the Killing form is negative definite (semisimple algebra).
    pub fn is_negative_definite(&self, n: usize) -> bool {
        let basis = SoAlgebra::basis(n);
        let b = self.matrix(&basis);
        // Check that all eigenvalues are negative
        let eigenvalues = b.symmetric_eigenvalues();
        for i in 0..eigenvalues.len() {
            if eigenvalues[i] >= 0.0 {
                return false;
            }
        }
        true
    }

    /// Check bilinearity.
    pub fn verify_bilinearity(&self, basis: &[LieAlgebraElement]) -> bool {
        if basis.len() < 2 {
            return true;
        }
        let x = &basis[0];
        let y = &basis[1];
        let a = 2.0;
        let b = 3.0;

        let ax = x.scale(a);
        let by = y.scale(b);
        let lhs = self.evaluate(&ax, &by, basis);
        let rhs = a * b * self.evaluate(x, y, basis);
        (lhs - rhs).abs() < 1e-8
    }

    /// Check symmetry: B(X,Y) = B(Y,X).
    pub fn verify_symmetry(&self, basis: &[LieAlgebraElement]) -> bool {
        let mat = self.matrix(basis);
        (mat - mat.transpose()).norm() < 1e-8
    }

    /// Check invariance: B([Z,X],Y) + B(X,[Z,Y]) = 0.
    pub fn verify_invariance(&self, basis: &[LieAlgebraElement]) -> bool {
        if basis.len() < 3 {
            return true;
        }
        let x = &basis[0];
        let y = &basis[1];
        let z = &basis[2];

        let zx = LieAlgebraElement::bracket(z, x);
        let zy = LieAlgebraElement::bracket(z, y);

        let term1 = self.evaluate(&zx, y, basis);
        let term2 = self.evaluate(x, &zy, basis);

        (term1 + term2).abs() < 1e-8
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_killing_so3() {
        let kf = KillingForm::for_so(3);
        assert_eq!(kf.algebra_dim, 3);
    }

    #[test]
    fn test_killing_so3_negative_definite() {
        let kf = KillingForm::for_so(3);
        assert!(kf.is_negative_definite(3));
    }

    #[test]
    fn test_killing_bilinearity() {
        let kf = KillingForm::for_so(3);
        let basis = SoAlgebra::basis(3);
        assert!(kf.verify_bilinearity(&basis));
    }

    #[test]
    fn test_killing_symmetry() {
        let kf = KillingForm::for_so(3);
        let basis = SoAlgebra::basis(3);
        assert!(kf.verify_symmetry(&basis));
    }

    #[test]
    fn test_killing_invariance() {
        let kf = KillingForm::for_so(3);
        let basis = SoAlgebra::basis(3);
        assert!(kf.verify_invariance(&basis));
    }

    #[test]
    fn test_killing_so2() {
        let kf = KillingForm::for_so(2);
        assert_eq!(kf.algebra_dim, 1);
    }

    #[test]
    fn test_killing_so4() {
        let kf = KillingForm::for_so(4);
        assert_eq!(kf.algebra_dim, 6);
    }

    #[test]
    fn test_killing_so4_negative_definite() {
        let kf = KillingForm::for_so(4);
        assert!(kf.is_negative_definite(4));
    }

    #[test]
    fn test_killing_so4_symmetry() {
        let kf = KillingForm::for_so(4);
        let basis = SoAlgebra::basis(4);
        assert!(kf.verify_symmetry(&basis));
    }

    #[test]
    fn test_killing_matrix() {
        let kf = KillingForm::for_so(3);
        let basis = SoAlgebra::basis(3);
        let mat = kf.matrix(&basis);
        assert_eq!(mat.nrows(), 3);
        assert_eq!(mat.ncols(), 3);
    }

    #[test]
    fn test_killing_so3_diagonal() {
        let kf = KillingForm::for_so(3);
        let basis = SoAlgebra::basis(3);
        let mat = kf.matrix(&basis);
        // For so(3), B is negative definite and proportional to identity
        // (up to normalization of basis)
        for i in 0..3 {
            assert!(mat[(i, i)] < 0.0);
        }
    }

    #[test]
    fn test_serde() {
        let kf = KillingForm::for_so(3);
        let json = serde_json::to_string(&kf).unwrap();
        let back: KillingForm = serde_json::from_str(&json).unwrap();
        assert_eq!(back.algebra_dim, 3);
    }
}
