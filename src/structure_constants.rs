//! Structure constants: [eᵢ, eⱼ] = Σ cᵢⱼᵏ eₖ.

use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};

use crate::lie_algebra::{LieAlgebraElement, SoAlgebra};

/// Structure constants tensor c^{k}_{ij} for a Lie algebra.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StructureConstants {
    /// Dimension of the Lie algebra.
    pub algebra_dim: usize,
    /// c[i][j][k] = c^k_{ij}
    pub constants: Vec<Vec<Vec<f64>>>,
}

impl StructureConstants {
    /// Compute structure constants for a given basis.
    pub fn from_basis(basis: &[LieAlgebraElement]) -> Self {
        let dim = basis.len();
        let n = basis[0].dim();

        let mut constants = vec![vec![vec![0.0; dim]; dim]; dim];

        for i in 0..dim {
            for j in 0..dim {
                let bracket = LieAlgebraElement::bracket(&basis[i], &basis[j]);
                // Decompose bracket into basis components
                for k in 0..dim {
                    let dot = (&bracket.matrix.component_mul(&basis[k].matrix)).sum();
                    let norm_sq = (&basis[k].matrix.component_mul(&basis[k].matrix)).sum();
                    if norm_sq.abs() > 1e-10 {
                        constants[i][j][k] = dot / norm_sq;
                    }
                }
            }
        }

        Self { algebra_dim: dim, constants }
    }

    /// Compute for so(n).
    pub fn for_so(n: usize) -> Self {
        let basis = SoAlgebra::basis(n);
        Self::from_basis(&basis)
    }

    /// Get c^k_{ij}.
    pub fn get(&self, i: usize, j: usize, k: usize) -> f64 {
        self.constants[i][j][k]
    }

    /// Verify antisymmetry: c^k_{ij} = -c^k_{ji}.
    pub fn verify_antisymmetry(&self) -> bool {
        for i in 0..self.algebra_dim {
            for j in 0..self.algebra_dim {
                for k in 0..self.algebra_dim {
                    if (self.constants[i][j][k] + self.constants[j][i][k]).abs() > 1e-8 {
                        return false;
                    }
                }
            }
        }
        true
    }

    /// Verify Jacobi identity from structure constants:
    /// c^m_{ij} c^l_{mk} + c^m_{jk} c^l_{mi} + c^m_{ki} c^l_{mj} = 0
    pub fn verify_jacobi(&self) -> bool {
        let d = self.algebra_dim;
        for i in 0..d {
            for j in 0..d {
                for k in 0..d {
                    for l in 0..d {
                        let mut sum = 0.0;
                        for m in 0..d {
                            sum += self.constants[i][j][m] * self.constants[m][k][l];
                            sum += self.constants[j][k][m] * self.constants[m][i][l];
                            sum += self.constants[k][i][m] * self.constants[m][j][l];
                        }
                        if sum.abs() > 1e-8 {
                            return false;
                        }
                    }
                }
            }
        }
        true
    }

    /// Dimension of the algebra.
    pub fn dim(&self) -> usize {
        self.algebra_dim
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_so3_structure_constants() {
        let sc = StructureConstants::for_so(3);
        assert_eq!(sc.dim(), 3);
    }

    #[test]
    fn test_so3_antisymmetry() {
        let sc = StructureConstants::for_so(3);
        assert!(sc.verify_antisymmetry());
    }

    #[test]
    fn test_so3_jacobi() {
        let sc = StructureConstants::for_so(3);
        assert!(sc.verify_jacobi());
    }

    #[test]
    fn test_so2_structure_constants() {
        let sc = StructureConstants::for_so(2);
        assert_eq!(sc.dim(), 1);
        // so(2) is abelian, all constants should be 0
        assert!(sc.constants[0][0][0].abs() < 1e-10);
    }

    #[test]
    fn test_so4_structure_constants() {
        let sc = StructureConstants::for_so(4);
        assert_eq!(sc.dim(), 6);
    }

    #[test]
    fn test_so4_jacobi() {
        let sc = StructureConstants::for_so(4);
        assert!(sc.verify_jacobi());
    }

    #[test]
    fn test_so4_antisymmetry() {
        let sc = StructureConstants::for_so(4);
        assert!(sc.verify_antisymmetry());
    }

    #[test]
    fn test_so3_levi_civita() {
        let sc = StructureConstants::for_so(3);
        // For so(3), c^k_{ij} = ε_{ijk} (Levi-Civita symbol)
        // [e1, e2] = e3 → c[0][1][2] = 1
        assert!((sc.get(0, 1, 2) - 1.0).abs() < 1e-8 ||
                (sc.get(0, 1, 2) + 1.0).abs() < 1e-8); // sign depends on normalization
    }

    #[test]
    fn test_so3_nonzero_bracket() {
        let sc = StructureConstants::for_so(3);
        // At least some constants should be nonzero
        let mut nonzero = 0;
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    if sc.get(i, j, k).abs() > 1e-8 {
                        nonzero += 1;
                    }
                }
            }
        }
        assert!(nonzero > 0);
    }

    #[test]
    fn test_from_custom_basis() {
        let e1 = LieAlgebraElement::new(DMatrix::from_vec(2, 2, vec![0.0, 1.0, -1.0, 0.0]));
        let sc = StructureConstants::from_basis(&[e1.clone()]);
        assert_eq!(sc.dim(), 1);
    }

    #[test]
    fn test_serde() {
        let sc = StructureConstants::for_so(3);
        let json = serde_json::to_string(&sc).unwrap();
        let back: StructureConstants = serde_json::from_str(&json).unwrap();
        assert_eq!(back.dim(), 3);
    }
}
