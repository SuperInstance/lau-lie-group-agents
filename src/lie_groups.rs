//! Classical Lie groups: SO(n), SU(n), GL(n), Sp(2n), SE(3).

use nalgebra::{DMatrix, DVector, ComplexField};
use serde::{Deserialize, Serialize};

/// A matrix Lie group element.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LieGroupElement {
    /// The matrix representation of the group element.
    pub matrix: DMatrix<f64>,
    /// Which group this belongs to.
    pub group_type: GroupType,
}

/// Types of Lie groups supported.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum GroupType {
    SO(usize),
    SU(usize),
    GL(usize),
    Sp(usize),
    SE3,
}

impl LieGroupElement {
    /// Create from a matrix.
    pub fn new(matrix: DMatrix<f64>, group_type: GroupType) -> Self {
        Self { matrix, group_type }
    }

    /// Identity element.
    pub fn identity(group_type: GroupType) -> Self {
        let dim = group_type.matrix_dim();
        Self {
            matrix: DMatrix::identity(dim, dim),
            group_type,
        }
    }

    /// Group multiplication.
    pub fn multiply(&self, other: &LieGroupElement) -> LieGroupElement {
        assert_eq!(self.group_type, other.group_type);
        LieGroupElement {
            matrix: &self.matrix * &other.matrix,
            group_type: self.group_type,
        }
    }

    /// Inverse.
    pub fn inverse(&self) -> Option<LieGroupElement> {
        let inv = self.matrix.clone().try_inverse()?;
        Some(LieGroupElement {
            matrix: inv,
            group_type: self.group_type,
        })
    }

    /// Matrix dimension.
    pub fn dim(&self) -> usize {
        self.matrix.nrows()
    }
}

impl GroupType {
    /// Matrix dimension for this group type.
    pub fn matrix_dim(&self) -> usize {
        match self {
            GroupType::SO(n) => *n,
            GroupType::SU(n) => *n,
            GroupType::GL(n) => *n,
            GroupType::Sp(n) => 2 * n,
            GroupType::SE3 => 4,
        }
    }
}

/// Special Orthogonal Group SO(n): n×n orthogonal matrices with det = 1.
pub struct SO;

impl SO {
    /// Identity in SO(n).
    pub fn identity(n: usize) -> LieGroupElement {
        LieGroupElement::identity(GroupType::SO(n))
    }

    /// Rotation in the (i,j)-plane by angle theta.
    pub fn rotation(n: usize, i: usize, j: usize, theta: f64) -> LieGroupElement {
        let mut m = DMatrix::identity(n, n);
        let c = theta.cos();
        let s = theta.sin();
        m[(i, i)] = c;
        m[(i, j)] = -s;
        m[(j, i)] = s;
        m[(j, j)] = c;
        LieGroupElement::new(m, GroupType::SO(n))
    }

    /// Check if a matrix is in SO(n).
    pub fn is_member(m: &DMatrix<f64>) -> bool {
        let n = m.nrows();
        if m.ncols() != n {
            return false;
        }
        let t_m = m.transpose();
        let product = &t_m * m;
        let id = DMatrix::identity(n, n);
        (product - id.clone()).norm() < 1e-8 && (m.determinant() - 1.0).abs() < 1e-8
    }

    /// SO(2) rotation by angle theta.
    pub fn so2(theta: f64) -> LieGroupElement {
        Self::rotation(2, 0, 1, theta)
    }

    /// SO(3) rotation around x-axis.
    pub fn rotx(theta: f64) -> LieGroupElement {
        Self::rotation(3, 1, 2, theta)
    }

    /// SO(3) rotation around y-axis.
    pub fn roty(theta: f64) -> LieGroupElement {
        Self::rotation(3, 0, 2, theta)
    }

    /// SO(3) rotation around z-axis.
    pub fn rotz(theta: f64) -> LieGroupElement {
        Self::rotation(3, 0, 1, theta)
    }
}

/// General Linear Group GL(n): all invertible n×n matrices.
pub struct GL;

impl GL {
    /// Identity in GL(n).
    pub fn identity(n: usize) -> LieGroupElement {
        LieGroupElement::identity(GroupType::GL(n))
    }

    /// Check if a matrix is in GL(n) (invertible).
    pub fn is_member(m: &DMatrix<f64>) -> bool {
        m.nrows() == m.ncols() && m.determinant().abs() > 1e-10
    }

    /// Create from a matrix.
    pub fn from_matrix(n: usize, m: DMatrix<f64>) -> LieGroupElement {
        LieGroupElement::new(m, GroupType::GL(n))
    }

    /// Diagonal matrix with given entries.
    pub fn diagonal(n: usize, entries: &[f64]) -> LieGroupElement {
        let mut m = DMatrix::zeros(n, n);
        for (i, &e) in entries.iter().enumerate().take(n) {
            m[(i, i)] = e;
        }
        LieGroupElement::new(m, GroupType::GL(n))
    }
}

/// Special Unitary Group SU(n): n×n unitary matrices with det = 1.
/// Represented as real 2n×2n matrices (realification).
pub struct SU;

impl SU {
    /// Identity in SU(n) as real matrix.
    pub fn identity(n: usize) -> LieGroupElement {
        LieGroupElement::identity(GroupType::SU(n))
    }

    /// U(1) ⊂ SU(2): rotation by angle theta, represented as 2×2 real.
    pub fn su2_rotation(theta: f64) -> LieGroupElement {
        let m = DMatrix::from_vec(2, 2, vec![
            theta.cos(), theta.sin(),
            -theta.sin(), theta.cos(),
        ]);
        LieGroupElement::new(m, GroupType::SU(2))
    }

    /// Check if a real 2n×2n matrix corresponds to an SU(n) element.
    pub fn is_member_real(m: &DMatrix<f64>) -> bool {
        let n = m.nrows();
        if m.ncols() != n || n % 2 != 0 {
            return false;
        }
        let t_m = m.transpose();
        let product = &t_m * m;
        let id = DMatrix::identity(n, n);
        (product - id.clone()).norm() < 1e-6 && (m.determinant() - 1.0).abs() < 1e-6
    }
}

/// Symplectic Group Sp(2n): preserves the symplectic form ω = J.
pub struct Sp;

impl Sp {
    /// Identity in Sp(2n).
    pub fn identity(n: usize) -> LieGroupElement {
        LieGroupElement::identity(GroupType::Sp(n))
    }

    /// The standard symplectic matrix J.
    pub fn j_matrix(n: usize) -> DMatrix<f64> {
        let dim = 2 * n;
        let mut j = DMatrix::zeros(dim, dim);
        for i in 0..n {
            j[(i, n + i)] = 1.0;
            j[(n + i, i)] = -1.0;
        }
        j
    }

    /// Check if a matrix is in Sp(2n): M^T J M = J.
    pub fn is_member(m: &DMatrix<f64>, n: usize) -> bool {
        let j = Self::j_matrix(n);
        let product = m.transpose() * &j * m;
        (product - j).norm() < 1e-6
    }

    /// A basic symplectic matrix: block diagonal with 2x2 rotations.
    pub fn symplectic_rotation(n: usize, theta: f64) -> LieGroupElement {
        let dim = 2 * n;
        let mut m = DMatrix::identity(dim, dim);
        let c = theta.cos();
        let s = theta.sin();
        // Each 2×2 block is a rotation
        for i in 0..n {
            m[(i, i)] = c;
            m[(i, n + i)] = -s;
            m[(n + i, i)] = s;
            m[(n + i, n + i)] = c;
        }
        LieGroupElement::new(m, GroupType::Sp(n))
    }
}

/// Special Euclidean Group SE(3): rigid body transformations.
pub struct SE3;

impl SE3 {
    /// Identity in SE(3).
    pub fn identity() -> LieGroupElement {
        LieGroupElement::identity(GroupType::SE3)
    }

    /// Create SE(3) element from rotation matrix and translation.
    pub fn from_rotation_translation(rotation: &DMatrix<f64>, translation: &DVector<f64>) -> LieGroupElement {
        assert_eq!(rotation.nrows(), 3);
        assert_eq!(rotation.ncols(), 3);
        assert_eq!(translation.len(), 3);
        let mut m = DMatrix::identity(4, 4);
        for i in 0..3 {
            for j in 0..3 {
                m[(i, j)] = rotation[(i, j)];
            }
            m[(i, 3)] = translation[i];
        }
        LieGroupElement::new(m, GroupType::SE3)
    }

    /// Pure translation.
    pub fn translation(tx: f64, ty: f64, tz: f64) -> LieGroupElement {
        let t = DVector::from_vec(vec![tx, ty, tz]);
        let r = DMatrix::identity(3, 3);
        Self::from_rotation_translation(&r, &t)
    }

    /// Extract rotation and translation.
    pub fn decompose(m: &DMatrix<f64>) -> (DMatrix<f64>, DVector<f64>) {
        let mut r = DMatrix::zeros(3, 3);
        let mut t = DVector::zeros(3);
        for i in 0..3 {
            for j in 0..3 {
                r[(i, j)] = m[(i, j)];
            }
            t[i] = m[(i, 3)];
        }
        (r, t)
    }

    /// Transform a 3D point.
    pub fn transform_point(m: &DMatrix<f64>, p: &DVector<f64>) -> DVector<f64> {
        let mut homo = DVector::zeros(4);
        homo[0] = p[0];
        homo[1] = p[1];
        homo[2] = p[2];
        homo[3] = 1.0;
        let result = m * &homo;
        DVector::from_vec(vec![result[0], result[1], result[2]])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_so_identity() {
        let id = SO::identity(3);
        assert_eq!(id.dim(), 3);
        assert!(SO::is_member(&id.matrix));
    }

    #[test]
    fn test_so2_rotation() {
        let r = SO::so2(std::f64::consts::PI / 2.0);
        assert!(SO::is_member(&r.matrix));
        assert!((r.matrix[(0, 0)] - 0.0).abs() < 1e-10);
        assert!((r.matrix[(0, 1)] - (-1.0)).abs() < 1e-10);
    }

    #[test]
    fn test_so3_rotations() {
        let rx = SO::rotx(std::f64::consts::PI);
        assert!(SO::is_member(&rx.matrix));
        let ry = SO::roty(std::f64::consts::PI);
        assert!(SO::is_member(&ry.matrix));
        let rz = SO::rotz(std::f64::consts::PI);
        assert!(SO::is_member(&rz.matrix));
    }

    #[test]
    fn test_so_inverse() {
        let r = SO::rotz(1.0);
        let inv = r.inverse().unwrap();
        let product = r.multiply(&inv);
        let id = SO::identity(3);
        assert!((product.matrix - id.matrix).norm() < 1e-10);
    }

    #[test]
    fn test_so_composition() {
        let r1 = SO::so2(1.0);
        let r2 = SO::so2(2.0);
        let r3 = r1.multiply(&r2);
        let r_direct = SO::so2(3.0);
        assert!((r3.matrix - r_direct.matrix).norm() < 1e-10);
    }

    #[test]
    fn test_gl_identity() {
        let id = GL::identity(3);
        assert!(GL::is_member(&id.matrix));
    }

    #[test]
    fn test_gl_diagonal() {
        let d = GL::diagonal(3, &[2.0, 3.0, 4.0]);
        assert!(GL::is_member(&d.matrix));
        assert!((d.matrix.determinant() - 24.0).abs() < 1e-10);
    }

    #[test]
    fn test_gl_not_member() {
        let m = DMatrix::zeros(2, 2);
        assert!(!GL::is_member(&m)); // det = 0
    }

    #[test]
    fn test_su2_rotation() {
        let r = SU::su2_rotation(1.0);
        assert!(SU::is_member_real(&r.matrix));
    }

    #[test]
    fn test_sp_identity() {
        let id = Sp::identity(2);
        assert!(Sp::is_member(&id.matrix, 2));
    }

    #[test]
    fn test_sp_j_matrix() {
        let j = Sp::j_matrix(1);
        assert!((j[(0, 1)] - 1.0).abs() < 1e-10);
        assert!((j[(1, 0)] - (-1.0)).abs() < 1e-10);
    }

    #[test]
    fn test_sp_rotation() {
        let r = Sp::symplectic_rotation(1, 0.5);
        assert!(Sp::is_member(&r.matrix, 1));
    }

    #[test]
    fn test_se3_identity() {
        let id = SE3::identity();
        assert_eq!(id.dim(), 4);
    }

    #[test]
    fn test_se3_translation() {
        let t = SE3::translation(1.0, 2.0, 3.0);
        let (r, trans) = SE3::decompose(&t.matrix);
        assert!((r - DMatrix::identity(3, 3)).norm() < 1e-10);
        assert!((trans[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_se3_transform_point() {
        let t = SE3::translation(1.0, 0.0, 0.0);
        let p = DVector::from_vec(vec![0.0, 0.0, 0.0]);
        let result = SE3::transform_point(&t.matrix, &p);
        assert!((result[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_se3_compose() {
        let t1 = SE3::translation(1.0, 0.0, 0.0);
        let t2 = SE3::translation(0.0, 2.0, 0.0);
        let t3 = t1.multiply(&t2);
        let (_, trans) = SE3::decompose(&t3.matrix);
        assert!((trans[0] - 1.0).abs() < 1e-10);
        assert!((trans[1] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_group_type_dim() {
        assert_eq!(GroupType::SO(3).matrix_dim(), 3);
        assert_eq!(GroupType::SU(2).matrix_dim(), 2);
        assert_eq!(GroupType::GL(4).matrix_dim(), 4);
        assert_eq!(GroupType::Sp(2).matrix_dim(), 4);
        assert_eq!(GroupType::SE3.matrix_dim(), 4);
    }

    #[test]
    fn test_serde() {
        let r = SO::rotz(1.0);
        let json = serde_json::to_string(&r).unwrap();
        let back: LieGroupElement = serde_json::from_str(&json).unwrap();
        assert_eq!(back.group_type, GroupType::SO(3));
    }
}
