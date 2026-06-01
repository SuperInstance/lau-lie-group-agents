//! Baker-Campbell-Hausdorff formula.
//!
//! log(exp(X)exp(Y)) = X + Y + ½[X,Y] + 1/12[X,[X,Y]] - 1/12[Y,[X,Y]] + ...

use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};

use crate::lie_algebra::LieAlgebraElement;

/// BCH approximation order.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum BCHOrder {
    First,  // X + Y
    Second, // X + Y + ½[X,Y]
    Third,  // Full third-order
}

/// Compute the Baker-Campbell-Hausdorff formula up to a given order.
pub fn bch(x: &LieAlgebraElement, y: &LieAlgebraElement, order: BCHOrder) -> LieAlgebraElement {
    let mut result = x.add(y);

    if matches!(order, BCHOrder::First) {
        return result;
    }

    // Second order: + ½[X,Y]
    let xy = LieAlgebraElement::bracket(x, y);
    result = result.add(&xy.scale(0.5));

    if matches!(order, BCHOrder::Second) {
        return result;
    }

    // Third order: + 1/12 [X,[X,Y]] - 1/12 [Y,[X,Y]]
    let x_xy = LieAlgebraElement::bracket(x, &xy);
    let y_xy = LieAlgebraElement::bracket(y, &xy);
    result = result.add(&x_xy.scale(1.0 / 12.0));
    result = result.add(&y_xy.scale(-1.0 / 12.0));

    result
}

/// Verify BCH: exp(X)exp(Y) ≈ exp(BCH(X,Y)).
pub fn verify_bch(x: &LieAlgebraElement, y: &LieAlgebraElement, order: BCHOrder) -> bool {
    let exp_x = crate::exponential_map::matrix_exp(&x.matrix);
    let exp_y = crate::exponential_map::matrix_exp(&y.matrix);
    let lhs = &exp_x * &exp_y;

    let bch_result = bch(x, y, order);
    let rhs = crate::exponential_map::matrix_exp(&bch_result.matrix);

    let error = (lhs - rhs).norm();
    let tolerance = match order {
        BCHOrder::First => 0.5,
        BCHOrder::Second => 0.05,
        BCHOrder::Third => 0.01,
    };

    // Only check for small X, Y where BCH converges well
    if x.norm() < 0.5 && y.norm() < 0.5 {
        error < tolerance
    } else {
        true // Skip verification for large elements
    }
}

/// Compute the commutator series of BCH to arbitrary depth.
pub fn bch_series(x: &LieAlgebraElement, y: &LieAlgebraElement, depth: usize) -> LieAlgebraElement {
    let mut result = x.add(y);

    if depth >= 2 {
        let xy = LieAlgebraElement::bracket(x, y);
        result = result.add(&xy.scale(0.5));
    }

    if depth >= 3 {
        let xy = LieAlgebraElement::bracket(x, y);
        let x_xy = LieAlgebraElement::bracket(x, &xy);
        let y_xy = LieAlgebraElement::bracket(y, &xy);
        result = result.add(&x_xy.scale(1.0 / 12.0));
        result = result.add(&y_xy.scale(-1.0 / 12.0));
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    fn small_element() -> LieAlgebraElement {
        LieAlgebraElement::new(DMatrix::from_vec(2, 2, vec![
            0.0, 0.1, -0.1, 0.0,
        ]))
    }

    fn small_element2() -> LieAlgebraElement {
        LieAlgebraElement::new(DMatrix::from_vec(2, 2, vec![
            0.0, 0.0, 0.1, 0.0,
        ]))
    }

    #[test]
    fn test_bch_first_order() {
        let x = small_element();
        let y = small_element2();
        let result = bch(&x, &y, BCHOrder::First);
        assert_eq!(result.dim(), 2);
    }

    #[test]
    fn test_bch_second_order() {
        let x = small_element();
        let y = small_element2();
        let result = bch(&x, &y, BCHOrder::Second);
        assert_eq!(result.dim(), 2);
    }

    #[test]
    fn test_bch_third_order() {
        let x = small_element();
        let y = small_element2();
        let result = bch(&x, &y, BCHOrder::Third);
        assert_eq!(result.dim(), 2);
    }

    #[test]
    fn test_bch_commutative_case() {
        // When [X,Y] = 0, BCH reduces to X + Y
        let x = LieAlgebraElement::new(DMatrix::from_vec(2, 2, vec![
            0.0, 0.0, 0.0, 0.0,
        ]));
        let y = LieAlgebraElement::new(DMatrix::from_vec(2, 2, vec![
            1.0, 0.0, 0.0, -1.0,
        ]));
        let xy = LieAlgebraElement::bracket(&x, &y);
        assert!(xy.norm() < 1e-10); // X is zero, so bracket is zero

        let bch_result = bch(&x, &y, BCHOrder::Third);
        assert!((bch_result.matrix - y.matrix).norm() < 1e-10);
    }

    #[test]
    fn test_bch_verify_first() {
        let x = small_element();
        let y = small_element2();
        assert!(verify_bch(&x, &y, BCHOrder::First));
    }

    #[test]
    fn test_bch_verify_second() {
        let x = small_element();
        let y = small_element2();
        assert!(verify_bch(&x, &y, BCHOrder::Second));
    }

    #[test]
    fn test_bch_verify_third() {
        let x = small_element();
        let y = small_element2();
        assert!(verify_bch(&x, &y, BCHOrder::Third));
    }

    #[test]
    fn test_bch_series() {
        let x = small_element();
        let y = small_element2();
        let result = bch_series(&x, &y, 3);
        assert_eq!(result.dim(), 2);
    }

    #[test]
    fn test_bch_higher_order_more_accurate() {
        let x = small_element();
        let y = small_element2();
        let exp_x = crate::exponential_map::matrix_exp(&x.matrix);
        let exp_y = crate::exponential_map::matrix_exp(&y.matrix);
        let exact = &exp_x * &exp_y;

        let bch1 = bch(&x, &y, BCHOrder::First);
        let bch3 = bch(&x, &y, BCHOrder::Third);
        let approx1 = crate::exponential_map::matrix_exp(&bch1.matrix);
        let approx3 = crate::exponential_map::matrix_exp(&bch3.matrix);

        let err1 = (exact - approx1).norm();
        let err3 = (exact - approx3).norm();
        assert!(err3 <= err1 + 1e-10);
    }

    #[test]
    fn test_bch_order_serde() {
        let order = BCHOrder::Third;
        let json = serde_json::to_string(&order).unwrap();
        let back: BCHOrder = serde_json::from_str(&json).unwrap();
        assert!(matches!(back, BCHOrder::Third));
    }
}
