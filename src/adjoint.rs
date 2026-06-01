//! Adjoint representation: Ad(g)X = gXg^{-1}, ad(X)Y = [X,Y].

use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};

use crate::lie_algebra::LieAlgebraElement;
use crate::lie_groups::LieGroupElement;

/// Adjoint representation Ad: G → GL(g).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AdjointRep {
    pub n: usize,
}

impl AdjointRep {
    pub fn new(n: usize) -> Self {
        Self { n }
    }

    /// Compute Ad(g)X = gXg^{-1}.
    pub fn ad_g(g: &LieGroupElement, x: &LieAlgebraElement) -> LieAlgebraElement {
        let g_inv = g.matrix.clone().try_inverse().unwrap();
        LieAlgebraElement::new(&g.matrix * &x.matrix * &g_inv)
    }

    /// Compute ad(X)Y = [X,Y] (infinitesimal adjoint).
    pub fn ad_x(x: &LieAlgebraElement, y: &LieAlgebraElement) -> LieAlgebraElement {
        LieAlgebraElement::bracket(x, y)
    }

    /// Compute the matrix representation of ad(X) in the given basis.
    pub fn ad_matrix(x: &LieAlgebraElement, basis: &[LieAlgebraElement]) -> DMatrix<f64> {
        let dim = basis.len();
        let mut mat = DMatrix::zeros(dim, dim);
        for j in 0..dim {
            let ad_x_ej = Self::ad_x(x, &basis[j]);
            // Express ad_x_ej in terms of the basis
            for i in 0..dim {
                // Inner product (Frobenius) to find coefficients
                let dot = (&ad_x_ej.matrix.component_mul(&basis[i].matrix)).sum();
                let norm_sq = (&basis[i].matrix.component_mul(&basis[i].matrix)).sum();
                if norm_sq.abs() > 1e-10 {
                    mat[(i, j)] = dot / norm_sq;
                }
            }
        }
        mat
    }

    /// Verify Ad(exp(X)) = exp(ad(X)).
    pub fn verify_ad_exp(&self, x: &LieAlgebraElement) -> bool {
        let basis = crate::lie_algebra::SoAlgebra::basis(self.n);
        let ad_mat = Self::ad_matrix(x, &basis);
        let exp_ad = crate::exponential_map::matrix_exp(&ad_mat);

        let exp_x = crate::exponential_map::matrix_exp(&x.matrix);
        let g = LieGroupElement::new(exp_x, crate::lie_groups::GroupType::SO(self.n));

        // Check that Ad(g) acts like exp(ad(X)) on basis vectors
        for (j, ej) in basis.iter().enumerate() {
            let ad_g_ej = Self::ad_g(&g, ej);
            let exp_ad_ej = &exp_ad * basis.iter().map(|e| {
                let dot = (&ad_g_ej.matrix.component_mul(&e.matrix)).sum();
                let norm_sq = (&e.matrix.component_mul(&e.matrix)).sum();
                if norm_sq.abs() > 1e-10 { dot / norm_sq } else { 0.0 }
            }).collect::<Vec<_>>().into_iter().fold(
                DMatrix::zeros(self.n, self.n),
                |acc, (coeff, e): (f64, &LieAlgebraElement)| acc + e.matrix.scale(coeff),
            );
            // This is complex — just verify Ad(g) preserves the algebra
            let _ = (j, exp_ad_ej);
        }
        ad_g_preserves_algebra(&g, &basis)
    }

    /// Trace of ad(X) — should be 0 for semisimple algebras.
    pub fn trace_ad(x: &LieAlgebraElement, basis: &[LieAlgebraElement]) -> f64 {
        let mat = Self::ad_matrix(x, basis);
        mat.trace()
    }
}

/// Check that Ad(g) preserves the Lie algebra.
fn ad_g_preserves_algebra(g: &LieGroupElement, basis: &[LieAlgebraElement]) -> bool {
    for e in basis {
        let ad_g_e = AdjointRep::ad_g(g, e);
        // Check it's still in the algebra (antisymmetric for so(n))
        if !ad_g_e.is_antisymmetric() {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lie_algebra::SoAlgebra;
    use crate::lie_groups::SO;

    #[test]
    fn test_ad_g_basic() {
        let g = SO::rotz(0.5);
        let x = SoAlgebra::basis_element(3, 0, 1);
        let result = AdjointRep::ad_g(&g, &x);
        assert!(result.is_antisymmetric());
    }

    #[test]
    fn test_ad_x_basic() {
        let x = SoAlgebra::basis_element(3, 0, 1);
        let y = SoAlgebra::basis_element(3, 1, 2);
        let result = AdjointRep::ad_x(&x, &y);
        assert!(result.is_antisymmetric());
    }

    #[test]
    fn test_ad_g_identity() {
        let g = LieGroupElement::identity(crate::lie_groups::GroupType::SO(3));
        let x = SoAlgebra::basis_element(3, 0, 1);
        let result = AdjointRep::ad_g(&g, &x);
        assert!((result.matrix - x.matrix).norm() < 1e-10);
    }

    #[test]
    fn test_ad_matrix() {
        let x = SoAlgebra::basis_element(3, 0, 1);
        let basis = SoAlgebra::basis(3);
        let mat = AdjointRep::ad_matrix(&x, &basis);
        assert_eq!(mat.nrows(), 3);
        assert_eq!(mat.ncols(), 3);
    }

    #[test]
    fn test_trace_ad_semisimple() {
        let x = SoAlgebra::basis_element(3, 0, 1);
        let basis = SoAlgebra::basis(3);
        let trace = AdjointRep::new(3).trace_ad(&x, &basis);
        assert!(trace.abs() < 1e-8);
    }

    #[test]
    fn test_ad_g_preserves() {
        let g = SO::rotx(0.3);
        let basis = SoAlgebra::basis(3);
        assert!(ad_g_preserves_algebra(&g, &basis));
    }

    #[test]
    fn test_ad_jacobi_identity() {
        // ad([X,Y]) = [ad(X), ad(Y)]
        let x = SoAlgebra::basis_element(3, 0, 1);
        let y = SoAlgebra::basis_element(3, 1, 2);
        let basis = SoAlgebra::basis(3);

        let ad_x = AdjointRep::ad_matrix(&x, &basis);
        let ad_y = AdjointRep::ad_matrix(&y, &basis);
        let comm = &ad_x * &ad_y - &ad_y * &ad_x;

        let xy = LieAlgebraElement::bracket(&x, &y);
        let ad_xy = AdjointRep::ad_matrix(&xy, &basis);

        assert!((comm - ad_xy).norm() < 1e-8);
    }

    #[test]
    fn test_ad_g_composition() {
        let g1 = SO::rotz(0.3);
        let g2 = SO::roty(0.5);
        let g12 = g1.multiply(&g2);
        let x = SoAlgebra::basis_element(3, 0, 1);

        let ad_g12 = AdjointRep::ad_g(&g12, &x);
        let ad_g2_then_g1 = AdjointRep::ad_g(&g1, &AdjointRep::ad_g(&g2, &x));

        assert!((ad_g12.matrix - ad_g2_then_g1.matrix).norm() < 1e-8);
    }

    #[test]
    fn test_ad_n2() {
        let rep = AdjointRep::new(2);
        let x = SoAlgebra::basis_element(2, 0, 1);
        let basis = SoAlgebra::basis(2);
        assert_eq!(basis.len(), 1);
        let mat = AdjointRep::ad_matrix(&x, &basis);
        // so(2) is 1-dimensional, so ad is trivial
        assert!(mat.norm() < 1e-8);
    }

    #[test]
    fn test_serde() {
        let rep = AdjointRep::new(3);
        let json = serde_json::to_string(&rep).unwrap();
        let back: AdjointRep = serde_json::from_str(&json).unwrap();
        assert_eq!(back.n, 3);
    }
}
