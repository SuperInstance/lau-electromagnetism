
/// Vacuum permittivity (F/m)
pub const EPSILON_0: f64 = 8.854187817e-12;

/// Vacuum permeability (H/m) = T·m/A
pub const MU_0: f64 = 4.0 * std::f64::consts::PI * 1e-7;

/// Speed of light in vacuum (m/s)
pub const C: f64 = 299_792_458.0;

/// Coulomb's constant (N·m²/C²)
pub const K_COULOMB: f64 = 1.0 / (4.0 * std::f64::consts::PI * EPSILON_0);

/// Elementary charge (C)
pub const E_CHARGE: f64 = 1.602176634e-19;

/// Electron mass (kg)
pub const ELECTRON_MASS: f64 = 9.1093837015e-31;

/// Proton mass (kg)
pub const PROTON_MASS: f64 = 1.67262192369e-27;

/// Boltzmann constant (J/K)
pub const BOLTZMANN: f64 = 1.380649e-23;

/// Characteristic impedance of free space (Ω)
pub const Z_0: f64 = MU_0 * C; // ≈ 376.73 Ω

/// Verify the fundamental relation: c = 1/sqrt(ε₀μ₀)
pub fn verify_speed_of_light() -> bool {
    let computed = 1.0 / (EPSILON_0 * MU_0).sqrt();
    (computed - C).abs() / C < 1e-6
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_speed_of_light_relation() {
        assert!(verify_speed_of_light());
    }

    #[test]
    fn test_coulomb_constant() {
        let expected = 8.987_551_787_368_176_4e9;
        assert!((K_COULOMB - expected).abs() / expected < 1e-6);
    }

    #[test]
    fn test_characteristic_impedance() {
        // Z₀ ≈ 376.73 Ω
        assert!((Z_0 - 376.73).abs() < 0.1);
    }
}
