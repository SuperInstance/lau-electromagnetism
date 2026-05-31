use crate::constants::{EPSILON_0, MU_0};

/// Material electromagnetic properties.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Material {
    pub name: String,
    pub relative_permittivity: f64,
    pub relative_permeability: f64,
    pub conductivity: f64,
}

impl Material {
    pub fn new(name: &str, epsilon_r: f64, mu_r: f64, sigma: f64) -> Self {
        Self { name: name.to_string(), relative_permittivity: epsilon_r, relative_permeability: mu_r, conductivity: sigma }
    }

    pub fn permittivity(&self) -> f64 { self.relative_permittivity * EPSILON_0 }
    pub fn permeability(&self) -> f64 { self.relative_permeability * MU_0 }
    pub fn wave_speed(&self) -> f64 { 1.0 / (self.permittivity() * self.permeability()).sqrt() }
    pub fn refractive_index(&self) -> f64 { (self.relative_permittivity * self.relative_permeability).sqrt() }
    pub fn intrinsic_impedance(&self) -> f64 { (self.permeability() / self.permittivity()).sqrt() }
    pub fn wavelength_in_material(&self, lambda0: f64) -> f64 { lambda0 / self.refractive_index() }

    pub fn skin_depth(&self, frequency: f64) -> f64 {
        if self.conductivity.abs() < f64::EPSILON || frequency.abs() < f64::EPSILON {
            return f64::INFINITY;
        }
        let omega = 2.0 * std::f64::consts::PI * frequency;
        (2.0 / (omega * self.permeability() * self.conductivity)).sqrt()
    }

    pub fn is_dielectric(&self) -> bool { self.conductivity < 1e-10 }
    pub fn is_conductor(&self) -> bool { self.conductivity > 1e3 }

    pub fn vacuum() -> Self { Self::new("vacuum", 1.0, 1.0, 0.0) }
    pub fn air() -> Self { Self::new("air", 1.0006, 1.0000004, 0.0) }
    pub fn glass() -> Self { Self::new("glass", 4.5, 1.0, 1e-12) }
    pub fn water() -> Self { Self::new("water", 80.0, 1.0, 5.5e-6) }
    pub fn copper() -> Self { Self::new("copper", 1.0, 0.999994, 5.96e7) }
    pub fn iron() -> Self { Self::new("iron", 1.0, 5000.0, 1.0e7) }
    pub fn teflon() -> Self { Self::new("teflon", 2.1, 1.0, 0.0) }
}

/// Electric displacement field: D = εE
pub fn displacement_field(epsilon: f64, e: f64) -> f64 { epsilon * e }

/// Magnetization: M = (μᵣ - 1)H
pub fn magnetization(mu_r: f64, h: f64) -> f64 { (mu_r - 1.0) * h }

/// B from H: B = μ₀μᵣH
pub fn b_from_h(mu_r: f64, h: f64) -> f64 { MU_0 * mu_r * h }

/// Polarization: P = ε₀(εᵣ - 1)E
pub fn polarization(epsilon_r: f64, e: f64) -> f64 { EPSILON_0 * (epsilon_r - 1.0) * e }

/// Tangential E continuity: E₁ₜ = E₂ₜ
pub fn tangential_e_continuity(e1_t: f64, e2_t: f64) -> bool { (e1_t - e2_t).abs() < 1e-10 }

/// Normal D continuity: ε₁E₁ₙ = ε₂E₂ₙ
pub fn normal_d_continuity(eps1: f64, e1_n: f64, eps2: f64, e2_n: f64) -> bool {
    (eps1 * e1_n - eps2 * e2_n).abs() < 1e-10
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::C;

    #[test]
    fn test_vacuum_properties() {
        let v = Material::vacuum();
        assert!((v.permittivity() - EPSILON_0).abs() < 1e-20);
        assert!((v.permeability() - MU_0).abs() < 1e-20);
        assert!((v.wave_speed() - C).abs() / C < 1e-6);
        assert!((v.refractive_index() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_glass_refractive_index() {
        let g = Material::glass();
        assert!((g.refractive_index() - 4.5_f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_water_high_permittivity() {
        assert!((Material::water().relative_permittivity - 80.0).abs() < 0.1);
    }

    #[test]
    fn test_copper_skin_depth() {
        let delta = Material::copper().skin_depth(1e6);
        assert!(delta > 0.0 && delta < 0.001);
    }

    #[test]
    fn test_iron_high_permeability() {
        assert!(Material::iron().relative_permeability > 1000.0);
    }

    #[test]
    fn test_teflon_low_loss() {
        let t = Material::teflon();
        assert!(t.is_dielectric());
        assert!(!t.is_conductor());
    }

    #[test]
    fn test_intrinsic_impedance_vacuum() {
        let eta = Material::vacuum().intrinsic_impedance();
        assert!((eta - 376.73).abs() < 0.1);
    }

    #[test]
    fn test_wavelength_in_glass() {
        let g = Material::glass();
        let lambda0 = 500e-9;
        assert!((g.wavelength_in_material(lambda0) - lambda0 / g.refractive_index()).abs() < 1e-20);
    }

    #[test]
    fn test_displacement_field() {
        assert!((displacement_field(EPSILON_0, 100.0) - EPSILON_0 * 100.0).abs() < 1e-15);
    }

    #[test]
    fn test_polarization() {
        let expected = EPSILON_0 * 3.0 * 100.0;
        assert!((polarization(4.0, 100.0) - expected).abs() < 1e-15);
    }
}
