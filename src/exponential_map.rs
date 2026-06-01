//! Matrix exponential map: exp: g → G.

use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};

use crate::lie_algebra::LieAlgebraElement;
use crate::lie_groups::LieGroupElement;

/// Compute the matrix exponential via scaling and squaring with Padé approximation.
pub fn matrix_exp(m: &DMatrix<f64>) -> DMatrix<f64> {
    let n = m.nrows();
    assert_eq!(n, m.ncols());

    // Use scaling and squaring with Taylor series
    let norm = m.norm();
    let scaling = if norm > 1.0 {
        (norm.ceil() as usize).next_power_of_two()
    } else {
        1
    };
    let scaling_f = scaling as f64;

    let scaled = m.scale(1.0 / scaling_f);

    // Taylor series: exp(A) ≈ I + A + A²/2! + ... + A^k/k!
    let mut result = DMatrix::identity(n, n);
    let mut term = DMatrix::identity(n, n);
    for k in 1..=30 {
        term = term * &scaled;
        for j in 0..n {
            for i in 0..n {
                term[(i, j)] /= k as f64;
            }
        }
        result += &term;
        if term.norm() < 1e-15 {
            break;
        }
    }

    // Squaring phase
    for _ in 0..scaling.trailing_zeros() as usize {
        result = &result * &result;
    }

    result
}

/// Compute the matrix logarithm (inverse of exp) via inverse scaling and squaring.
pub fn matrix_log(m: &DMatrix<f64>) -> DMatrix<f64> {
    let n = m.nrows();
    assert_eq!(n, m.ncols());

    // Use the iterative method: log(M) = (M-I) - (M-I)²/2 + ...
    // First, scale M close to identity
    let mut scaled = m.clone();
    let mut num_squarings = 0;
    while (&scaled - DMatrix::<f64>::identity(n, n)).norm() > 0.5 {
        scaled = scaled.map(|x| x.sqrt());
        num_squarings += 1;
    }

    let diff = &scaled - DMatrix::<f64>::identity(n, n);

    // Taylor series for log(I + X)
    let mut result = DMatrix::zeros(n, n);
    let mut term = diff.clone();
    for k in 1..=50 {
        let sign = if k % 2 == 0 { -1.0 } else { 1.0 };
        result += term.scale(sign / k as f64);
        term = &term * &diff;
        if term.norm() < 1e-15 {
            break;
        }
    }

    // Undo squarings: log(M) = 2^s * log(M^{1/2^s})
    result.scale((1 << num_squarings) as f64)
}

/// The exponential map from a Lie algebra to the Lie group.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExponentialMap {
    pub n: usize,
}

impl ExponentialMap {
    pub fn new(n: usize) -> Self {
        Self { n }
    }

    /// Map an algebra element to the group: exp(X).
    pub fn apply(&self, x: &LieAlgebraElement) -> LieGroupElement {
        let exp_x = matrix_exp(&x.matrix);
        LieGroupElement::new(exp_x, crate::lie_groups::GroupType::GL(self.n))
    }

    /// Verify exp(0) = I.
    pub fn verify_identity(&self) -> bool {
        let zero = LieAlgebraElement::zero(self.n);
        let result = self.apply(&zero);
        let id = DMatrix::identity(self.n, self.n);
        (result.matrix - id).norm() < 1e-10
    }

    /// Verify exp maps so(n) to SO(n).
    pub fn verify_so_mapping(&self) -> bool {
        let basis = crate::lie_algebra::SoAlgebra::basis(self.n);
        for x in &basis {
            let g = self.apply(x);
            if !crate::lie_groups::SO::is_member(&g.matrix) {
                return false;
            }
        }
        true
    }

    /// Compute exp(tX) for scalar t.
    pub fn exp_t(&self, x: &LieAlgebraElement, t: f64) -> DMatrix<f64> {
        matrix_exp(&x.matrix.scale(t))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lie_algebra::SoAlgebra;

    #[test]
    fn test_matrix_exp_zero() {
        let zero = DMatrix::zeros(3, 3);
        let result = matrix_exp(&zero);
        assert!((result - DMatrix::identity(3, 3)).norm() < 1e-10);
    }

    #[test]
    fn test_matrix_exp_identity() {
        let id = DMatrix::identity(2, 2);
        let result = matrix_exp(&id);
        let expected = DMatrix::from_vec(2, 2, vec![
            std::f64::consts::E, 0.0, 0.0, std::f64::consts::E,
        ]);
        assert!((result - expected).norm() < 1e-8);
    }

    #[test]
    fn test_matrix_exp_antisymmetric() {
        let x = SoAlgebra::basis_element(2, 0, 1);
        let result = matrix_exp(&x.matrix);
        assert!(crate::lie_groups::SO::is_member(&result));
    }

    #[test]
    fn test_matrix_exp_so3() {
        let x = SoAlgebra::basis_element(3, 0, 1);
        let result = matrix_exp(&x.matrix);
        assert!(crate::lie_groups::SO::is_member(&result));
    }

    #[test]
    fn test_exp_map_identity() {
        let map = ExponentialMap::new(3);
        assert!(map.verify_identity());
    }

    #[test]
    fn test_exp_map_so() {
        let map = ExponentialMap::new(3);
        assert!(map.verify_so_mapping());
    }

    #[test]
    fn test_exp_t() {
        let x = SoAlgebra::basis_element(2, 0, 1);
        let map = ExponentialMap::new(2);
        let result = map.exp_t(&x, 0.0);
        assert!((result - DMatrix::identity(2, 2)).norm() < 1e-10);
    }

    #[test]
    fn test_exp_so2_full_rotation() {
        let x = SoAlgebra::basis_element(2, 0, 1);
        let result = matrix_exp(&x.matrix.scale(2.0 * std::f64::consts::PI));
        let id = DMatrix::identity(2, 2);
        assert!((result - id).norm() < 1e-8);
    }

    #[test]
    fn test_matrix_log_roundtrip() {
        let m = DMatrix::from_vec(2, 2, vec![2.0, 0.0, 0.0, 3.0]);
        let log_m = matrix_log(&m);
        let exp_log = matrix_exp(&log_m);
        assert!((exp_log - m).norm() < 1e-6);
    }

    #[test]
    fn test_serde() {
        let map = ExponentialMap::new(3);
        let json = serde_json::to_string(&map).unwrap();
        let back: ExponentialMap = serde_json::from_str(&json).unwrap();
        assert_eq!(back.n, 3);
    }
}
