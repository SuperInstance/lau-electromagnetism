use crate::constants::{EPSILON_0, MU_0};
use crate::vector::Vec3;

/// Maxwell's equations in differential form.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct MaxwellFields {
    pub e_field: Vec3,
    pub b_field: Vec3,
    pub j_current: Vec3,
    pub rho: f64,
    pub d_e_dt: Vec3,
    pub d_b_dt: Vec3,
}

impl MaxwellFields {
    /// Gauss's law for electricity: ∇·E = ρ/ε₀
    pub fn gauss_law_divergence(&self, div_e: f64) -> f64 {
        div_e - self.rho / EPSILON_0
    }

    /// Gauss's law for magnetism: ∇·B = 0
    pub fn gauss_magnetism(&self, div_b: f64) -> f64 {
        div_b
    }

    /// Faraday's law: ∇×E = -∂B/∂t
    pub fn faradays_law(&self, curl_e: Vec3) -> Vec3 {
        curl_e + self.d_b_dt
    }

    /// Ampere-Maxwell law: ∇×B = μ₀J + μ₀ε₀ ∂E/∂t
    pub fn ampere_maxwell_law(&self, curl_b: Vec3) -> Vec3 {
        let expected = MU_0 * self.j_current + MU_0 * EPSILON_0 * self.d_e_dt;
        curl_b - expected
    }
}

/// Displacement current density: J_d = ε₀ ∂E/∂t
pub fn displacement_current(d_e_dt: Vec3) -> Vec3 {
    EPSILON_0 * d_e_dt
}

/// Gauss's law (integral): ∮ E·dA = Q_enc / ε₀
pub fn gauss_law_integral(e_field_magnitude: f64, area: f64, q_enclosed: f64) -> f64 {
    e_field_magnitude * area - q_enclosed / EPSILON_0
}

/// Faraday's law (integral): ∮ E·dl = -dΦ_B/dt
pub fn faradays_law_integral(emf: f64, d_flux_dt: f64) -> f64 {
    emf + d_flux_dt
}

/// Ampere's law (integral): ∮ B·dl = μ₀ I_enc + μ₀ε₀ dΦ_E/dt
pub fn amperes_law_integral(circulation: f64, i_enc: f64, d_flux_e_dt: f64) -> f64 {
    circulation - MU_0 * i_enc - MU_0 * EPSILON_0 * d_flux_e_dt
}

/// Induced EMF: EMF = -dΦ_B/dt
pub fn induced_emf(d_flux_dt: f64) -> f64 {
    -d_flux_dt
}

/// Magnetic flux through an area: Φ = B·A
pub fn magnetic_flux(b: f64, area: f64) -> f64 {
    b * area
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::vector;

    #[test]
    fn test_gauss_law_integral() {
        let q = 1e-6;
        let r = 1.0;
        let e = crate::constants::K_COULOMB * q / (r * r);
        let area = 4.0 * std::f64::consts::PI * r * r;
        let residual = gauss_law_integral(e, area, q);
        assert!(residual.abs() / (q / EPSILON_0) < 1e-6);
    }

    #[test]
    fn test_faradays_law_emf() {
        let d_flux_dt = 0.5;
        let emf = induced_emf(d_flux_dt);
        assert!((emf - (-0.5)).abs() < 1e-10);
    }

    #[test]
    fn test_magnetic_flux() {
        let flux = magnetic_flux(0.5, 2.0);
        assert!((flux - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_displacement_current() {
        let d_e_dt = vector![1e10, 0.0, 0.0];
        let j_d = displacement_current(d_e_dt);
        let expected = EPSILON_0 * 1e10;
        assert!((j_d.x - expected).abs() / expected < 1e-6);
    }

    #[test]
    fn test_maxwell_fields_gauss_electric() {
        let fields = MaxwellFields {
            e_field: Vec3::zeros(),
            b_field: Vec3::zeros(),
            j_current: Vec3::zeros(),
            rho: 1.0,
            d_e_dt: Vec3::zeros(),
            d_b_dt: Vec3::zeros(),
        };
        let div_e = 1.0 / EPSILON_0;
        let residual = fields.gauss_law_divergence(div_e);
        assert!(residual.abs() < 1e-10);
    }

    #[test]
    fn test_faradays_law_differential() {
        let fields = MaxwellFields {
            e_field: Vec3::zeros(),
            b_field: Vec3::zeros(),
            j_current: Vec3::zeros(),
            rho: 0.0,
            d_e_dt: Vec3::zeros(),
            d_b_dt: vector![0.0, 0.0, -1.0],
        };
        let curl_e = vector![0.0, 0.0, 1.0];
        let residual = fields.faradays_law(curl_e);
        assert!(residual.norm() < 1e-10);
    }
}
