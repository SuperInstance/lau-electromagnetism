use crate::constants::{EPSILON_0, K_COULOMB};
use crate::vector::{distance, unit_vector, Vec3};

/// A point charge with position and magnitude.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PointCharge {
    pub position: Vec3,
    pub charge: f64,
}

impl PointCharge {
    pub fn new(position: Vec3, charge: f64) -> Self {
        Self { position, charge }
    }
}

/// Compute Coulomb force on charge q2 due to charge q1.
/// F = k * q1 * q2 / r² * r̂₁₂
pub fn coulomb_force(q1: &PointCharge, q2: &PointCharge) -> Vec3 {
    let r = distance(q1.position, q2.position);
    if r < f64::EPSILON {
        return Vec3::zeros();
    }
    let magnitude = K_COULOMB * q1.charge * q2.charge / (r * r);
    // Force direction: from q1 to q2 (repulsive for same-sign charges)
    let direction = unit_vector(q1.position, q2.position);
    magnitude * direction
}

/// Compute electric field at point `p` due to a point charge.
/// E = k * q / r² * r̂
pub fn electric_field_point(charge: &PointCharge, point: Vec3) -> Vec3 {
    let r = distance(charge.position, point);
    if r < f64::EPSILON {
        return Vec3::zeros();
    }
    let magnitude = K_COULOMB * charge.charge / (r * r);
    let direction = unit_vector(charge.position, point);
    magnitude * direction
}

/// Compute total electric field at a point from multiple charges (superposition).
pub fn electric_field_superposition(charges: &[PointCharge], point: Vec3) -> Vec3 {
    charges.iter().map(|c| electric_field_point(c, point)).sum()
}

/// Compute electric flux through a sphere of given radius centered at origin,
/// enclosing a set of charges. This is Gauss's law in integral form.
/// Φ = Q_enc / ε₀
pub fn electric_flux_sphere(charges: &[PointCharge], radius: f64) -> f64 {
    let enclosed: f64 = charges
        .iter()
        .filter(|c| distance(c.position, Vec3::zeros()) < radius)
        .map(|c| c.charge)
        .sum();
    enclosed / EPSILON_0
}

/// Compute enclosed charge within a sphere of given radius.
pub fn enclosed_charge(charges: &[PointCharge], radius: f64) -> f64 {
    charges
        .iter()
        .filter(|c| distance(c.position, Vec3::zeros()) < radius)
        .map(|c| c.charge)
        .sum()
}

/// Electric field from an infinite line charge at distance d.
/// E = λ / (2π ε₀ d)
pub fn electric_field_line_charge(lambda: f64, d: f64) -> f64 {
    lambda / (2.0 * std::f64::consts::PI * EPSILON_0 * d)
}

/// Electric field from an infinite plane charge.
/// E = σ / (2ε₀)
pub fn electric_field_plane_charge(sigma: f64) -> f64 {
    sigma / (2.0 * EPSILON_0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::vector;

    #[test]
    fn test_coulomb_law_positive_charges_repel() {
        let q1 = PointCharge::new(vector![0.0, 0.0, 0.0], 1e-6);
        let q2 = PointCharge::new(vector![1.0, 0.0, 0.0], 1e-6);
        let force = coulomb_force(&q1, &q2);
        assert!(force.x > 0.0);
    }

    #[test]
    fn test_coulomb_law_opposite_charges_attract() {
        let q1 = PointCharge::new(vector![0.0, 0.0, 0.0], 1e-6);
        let q2 = PointCharge::new(vector![1.0, 0.0, 0.0], -1e-6);
        let force = coulomb_force(&q1, &q2);
        assert!(force.x < 0.0);
    }

    #[test]
    fn test_coulomb_law_magnitude() {
        let q1 = PointCharge::new(vector![0.0, 0.0, 0.0], 1e-6);
        let q2 = PointCharge::new(vector![1.0, 0.0, 0.0], 1e-6);
        let force = coulomb_force(&q1, &q2);
        let expected = K_COULOMB * 1e-12;
        assert!((force.norm() - expected).abs() / expected < 1e-6);
    }

    #[test]
    fn test_electric_field_positive_charge() {
        let q = PointCharge::new(vector![0.0, 0.0, 0.0], 1e-6);
        let p = vector![1.0, 0.0, 0.0];
        let field = electric_field_point(&q, p);
        assert!(field.x > 0.0);
    }

    #[test]
    fn test_electric_field_negative_charge() {
        let q = PointCharge::new(vector![0.0, 0.0, 0.0], -1e-6);
        let p = vector![1.0, 0.0, 0.0];
        let field = electric_field_point(&q, p);
        assert!(field.x < 0.0);
    }

    #[test]
    fn test_superposition() {
        let q1 = PointCharge::new(vector![1.0, 0.0, 0.0], 1e-6);
        let q2 = PointCharge::new(vector![-1.0, 0.0, 0.0], 1e-6);
        let p = vector![0.0, 0.0, 0.0];
        let field = electric_field_superposition(&[q1, q2], p);
        assert!(field.x.abs() < 1e-10);
    }

    #[test]
    fn test_gauss_law_sphere() {
        let charges = vec![
            PointCharge::new(vector![0.0, 0.0, 0.0], 1e-6),
            PointCharge::new(vector![0.5, 0.0, 0.0], 2e-6),
            PointCharge::new(vector![5.0, 0.0, 0.0], 3e-6),
        ];
        let flux = electric_flux_sphere(&charges, 1.0);
        let expected = 3e-6 / EPSILON_0;
        assert!((flux - expected).abs() / expected < 1e-6);
    }

    #[test]
    fn test_line_charge_field() {
        let e = electric_field_line_charge(1e-6, 0.1);
        assert!(e > 0.0);
        let e2 = electric_field_line_charge(1e-6, 0.2);
        assert!((e / e2 - 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_plane_charge_field() {
        let e = electric_field_plane_charge(1e-6);
        let expected = 1e-6 / (2.0 * EPSILON_0);
        assert!((e - expected).abs() / expected < 1e-6);
    }

    #[test]
    fn test_coulomb_inverse_square() {
        let q = PointCharge::new(vector![0.0, 0.0, 0.0], 1.0);
        let e1 = electric_field_point(&q, vector![1.0, 0.0, 0.0]).norm();
        let e2 = electric_field_point(&q, vector![2.0, 0.0, 0.0]).norm();
        assert!((e1 / e2 - 4.0).abs() < 1e-6);
    }
}
