//! Root systems: Aₙ, Bₙ, Cₙ, Dₙ.

use nalgebra::DVector;
use serde::{Deserialize, Serialize};

/// Classical root system type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum RootSystemType {
    A, // su(n+1)
    B, // so(2n+1)
    C, // sp(2n)
    D, // so(2n)
}

/// A root system with roots and their properties.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RootSystem {
    pub system_type: RootSystemType,
    pub rank: usize,
    /// Positive roots as vectors in R^rank.
    pub positive_roots: Vec<DVector<f64>>,
    /// Simple roots.
    pub simple_roots: Vec<DVector<f64>>,
}

impl RootSystem {
    /// Construct the root system Aₙ (associated to su(n+1)).
    pub fn type_a(n: usize) -> Self {
        if n == 0 {
            return Self {
                system_type: RootSystemType::A,
                rank: 0,
                positive_roots: vec![],
                simple_roots: vec![],
            };
        }

        // Simple roots: α_i = e_i - e_{i+1} for i = 1,...,n
        let simple_roots: Vec<DVector<f64>> = (0..n)
            .map(|i| {
                let mut v = DVector::zeros(n);
                v[i] = 1.0;
                if i + 1 < n {
                    v[i + 1] = -1.0;
                }
                v
            })
            .collect();

        // Positive roots: e_i - e_j for i < j
        let mut positive_roots = Vec::new();
        for i in 0..n {
            for j in (i + 1)..=n {
                let mut v = DVector::zeros(n + 1);
                v[i] = 1.0;
                v[j] = -1.0;
                // Project to rank n
                positive_roots.push(v.rows(0, n).into());
            }
        }

        Self { system_type: RootSystemType::A, rank: n, positive_roots, simple_roots }
    }

    /// Construct the root system Bₙ (associated to so(2n+1)).
    pub fn type_b(n: usize) -> Self {
        if n == 0 {
            return Self {
                system_type: RootSystemType::B,
                rank: 0,
                positive_roots: vec![],
                simple_roots: vec![],
            };
        }

        let simple_roots: Vec<DVector<f64>> = (0..n).map(|i| {
            let mut v = DVector::zeros(n);
            v[i] = 1.0;
            if i + 1 < n {
                v[i + 1] = -1.0;
            } else {
                // Last simple root: e_n
                v[i] = 1.0;
            }
            v
        }).collect();

        let mut positive_roots = Vec::new();
        // e_i ± e_j (i < j)
        for i in 0..n {
            for j in (i + 1)..n {
                let mut v1 = DVector::zeros(n);
                v1[i] = 1.0;
                v1[j] = 1.0;
                positive_roots.push(v1);

                let mut v2 = DVector::zeros(n);
                v2[i] = 1.0;
                v2[j] = -1.0;
                positive_roots.push(v2);
            }
        }
        // e_i
        for i in 0..n {
            let mut v = DVector::zeros(n);
            v[i] = 1.0;
            positive_roots.push(v);
        }

        Self { system_type: RootSystemType::B, rank: n, positive_roots, simple_roots }
    }

    /// Construct the root system Cₙ (associated to sp(2n)).
    pub fn type_c(n: usize) -> Self {
        if n == 0 {
            return Self {
                system_type: RootSystemType::C,
                rank: 0,
                positive_roots: vec![],
                simple_roots: vec![],
            };
        }

        let simple_roots: Vec<DVector<f64>> = (0..n).map(|i| {
            let mut v = DVector::zeros(n);
            v[i] = 1.0;
            if i + 1 < n {
                v[i + 1] = -1.0;
            }
            v
        }).collect();

        let mut positive_roots = Vec::new();
        // e_i ± e_j (i < j)
        for i in 0..n {
            for j in (i + 1)..n {
                let mut v1 = DVector::zeros(n);
                v1[i] = 1.0;
                v1[j] = 1.0;
                positive_roots.push(v1);

                let mut v2 = DVector::zeros(n);
                v2[i] = 1.0;
                v2[j] = -1.0;
                positive_roots.push(v2);
            }
        }
        // 2e_i
        for i in 0..n {
            let mut v = DVector::zeros(n);
            v[i] = 2.0;
            positive_roots.push(v);
        }

        Self { system_type: RootSystemType::C, rank: n, positive_roots, simple_roots }
    }

    /// Construct the root system Dₙ (associated to so(2n)).
    pub fn type_d(n: usize) -> Self {
        if n < 2 {
            return Self {
                system_type: RootSystemType::D,
                rank: n,
                positive_roots: vec![],
                simple_roots: vec![],
            };
        }

        let simple_roots: Vec<DVector<f64>> = (0..n).map(|i| {
            let mut v = DVector::zeros(n);
            if i < n - 1 {
                v[i] = 1.0;
                if i + 1 < n - 1 {
                    v[i + 1] = -1.0;
                }
            } else {
                // Last simple root: e_{n-1} + e_n
                v[n - 2] = 1.0;
                v[n - 1] = 1.0;
            }
            v
        }).collect();

        let mut positive_roots = Vec::new();
        // e_i ± e_j (i < j)
        for i in 0..n {
            for j in (i + 1)..n {
                let mut v1 = DVector::zeros(n);
                v1[i] = 1.0;
                v1[j] = 1.0;
                positive_roots.push(v1);

                let mut v2 = DVector::zeros(n);
                v2[i] = 1.0;
                v2[j] = -1.0;
                positive_roots.push(v2);
            }
        }

        Self { system_type: RootSystemType::D, rank: n, positive_roots, simple_roots }
    }

    /// Number of positive roots.
    pub fn num_positive_roots(&self) -> usize {
        self.positive_roots.len()
    }

    /// Number of simple roots.
    pub fn num_simple_roots(&self) -> usize {
        self.simple_roots.len()
    }

    /// Total number of roots (positive + negative).
    pub fn num_roots(&self) -> usize {
        2 * self.positive_roots.len()
    }

    /// Compute the Cartan matrix from simple roots.
    pub fn cartan_matrix(&self) -> nalgebra::DMatrix<f64> {
        let n = self.simple_roots.len();
        let mut c = nalgebra::DMatrix::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                let dot_ij = self.simple_roots[i].dot(&self.simple_roots[j]);
                let dot_ii = self.simple_roots[i].dot(&self.simple_roots[i]);
                if dot_ii.abs() > 1e-10 {
                    c[(i, j)] = 2.0 * dot_ij / dot_ii;
                }
            }
        }
        c
    }

    /// Verify Cartan matrix has 2's on diagonal.
    pub fn verify_cartan_diagonal(&self) -> bool {
        let c = self.cartan_matrix();
        for i in 0..c.nrows() {
            if (c[(i, i)] - 2.0).abs() > 1e-8 {
                return false;
            }
        }
        true
    }
}

/// Count positive roots for classical types.
pub fn count_positive_roots(system_type: RootSystemType, rank: usize) -> usize {
    match system_type {
        RootSystemType::A => rank * (rank + 1) / 2,
        RootSystemType::B => rank * rank,
        RootSystemType::C => rank * rank,
        RootSystemType::D => rank * (rank - 1),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_a2_roots() {
        let rs = RootSystem::type_a(2);
        assert_eq!(rs.num_simple_roots(), 2);
        assert_eq!(rs.num_positive_roots(), 3);
    }

    #[test]
    fn test_a3_roots() {
        let rs = RootSystem::type_a(3);
        assert_eq!(rs.num_simple_roots(), 3);
        assert_eq!(rs.num_positive_roots(), 6);
    }

    #[test]
    fn test_b2_roots() {
        let rs = RootSystem::type_b(2);
        assert_eq!(rs.num_positive_roots(), 4);
    }

    #[test]
    fn test_c3_roots() {
        let rs = RootSystem::type_c(3);
        assert_eq!(rs.num_positive_roots(), 9);
    }

    #[test]
    fn test_d3_roots() {
        let rs = RootSystem::type_d(3);
        assert_eq!(rs.num_positive_roots(), 6);
    }

    #[test]
    fn test_d4_roots() {
        let rs = RootSystem::type_d(4);
        assert_eq!(rs.num_positive_roots(), 12);
    }

    #[test]
    fn test_a2_cartan_matrix() {
        let rs = RootSystem::type_a(2);
        let c = rs.cartan_matrix();
        assert_eq!(c.nrows(), 2);
        assert!((c[(0, 0)] - 2.0).abs() < 1e-8);
        assert!((c[(1, 1)] - 2.0).abs() < 1e-8);
    }

    #[test]
    fn test_cartan_diagonal() {
        let rs = RootSystem::type_a(3);
        assert!(rs.verify_cartan_diagonal());
    }

    #[test]
    fn test_count_a() {
        assert_eq!(count_positive_roots(RootSystemType::A, 3), 6);
    }

    #[test]
    fn test_count_b() {
        assert_eq!(count_positive_roots(RootSystemType::B, 3), 9);
    }

    #[test]
    fn test_count_c() {
        assert_eq!(count_positive_roots(RootSystemType::C, 2), 4);
    }

    #[test]
    fn test_count_d() {
        assert_eq!(count_positive_roots(RootSystemType::D, 4), 12);
    }

    #[test]
    fn test_total_roots() {
        let rs = RootSystem::type_a(2);
        assert_eq!(rs.num_roots(), 6); // 3 positive + 3 negative
    }

    #[test]
    fn test_a1_simple_root() {
        let rs = RootSystem::type_a(1);
        assert_eq!(rs.num_simple_roots(), 1);
        assert_eq!(rs.num_positive_roots(), 1);
    }

    #[test]
    fn test_serde() {
        let rs = RootSystem::type_a(2);
        let json = serde_json::to_string(&rs).unwrap();
        let back: RootSystem = serde_json::from_str(&json).unwrap();
        assert_eq!(back.rank, 2);
    }
}
