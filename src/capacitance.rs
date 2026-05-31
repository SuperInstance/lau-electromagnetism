use crate::constants::{EPSILON_0, MU_0};

/// Parallel plate capacitor.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ParallelPlateCapacitor {
    pub area: f64,
    pub separation: f64,
    pub dielectric_k: f64,
}

impl ParallelPlateCapacitor {
    pub fn new(area: f64, separation: f64) -> Self {
        Self { area, separation, dielectric_k: 1.0 }
    }

    pub fn with_dielectric(area: f64, separation: f64, k: f64) -> Self {
        Self { area, separation, dielectric_k: k }
    }

    pub fn capacitance(&self) -> f64 {
        EPSILON_0 * self.dielectric_k * self.area / self.separation
    }

    pub fn energy(&self, voltage: f64) -> f64 {
        0.5 * self.capacitance() * voltage * voltage
    }

    pub fn electric_field(&self, voltage: f64) -> f64 {
        voltage / self.separation
    }

    pub fn charge(&self, voltage: f64) -> f64 {
        self.capacitance() * voltage
    }
}

/// Spherical capacitor.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct SphericalCapacitor {
    pub inner_radius: f64,
    pub outer_radius: f64,
    pub dielectric_k: f64,
}

impl SphericalCapacitor {
    pub fn new(a: f64, b: f64) -> Self {
        Self { inner_radius: a, outer_radius: b, dielectric_k: 1.0 }
    }

    pub fn capacitance(&self) -> f64 {
        4.0 * std::f64::consts::PI * EPSILON_0 * self.dielectric_k
            * self.inner_radius * self.outer_radius
            / (self.outer_radius - self.inner_radius)
    }
}

/// Cylindrical capacitor.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct CylindricalCapacitor {
    pub inner_radius: f64,
    pub outer_radius: f64,
    pub length: f64,
    pub dielectric_k: f64,
}

impl CylindricalCapacitor {
    pub fn new(a: f64, b: f64, length: f64) -> Self {
        Self { inner_radius: a, outer_radius: b, length, dielectric_k: 1.0 }
    }

    pub fn capacitance(&self) -> f64 {
        2.0 * std::f64::consts::PI * EPSILON_0 * self.dielectric_k * self.length
            / (self.outer_radius / self.inner_radius).ln()
    }
}

/// Solenoid inductor.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct SolenoidInductor {
    pub turns: f64,
    pub length: f64,
    pub area: f64,
    pub core_mu_r: f64,
}

impl SolenoidInductor {
    pub fn new(turns: f64, length: f64, area: f64) -> Self {
        Self { turns, length, area, core_mu_r: 1.0 }
    }

    pub fn with_core(turns: f64, length: f64, area: f64, mu_r: f64) -> Self {
        Self { turns, length, area, core_mu_r: mu_r }
    }

    pub fn inductance(&self) -> f64 {
        MU_0 * self.core_mu_r * self.turns * self.turns * self.area / self.length
    }

    pub fn energy(&self, current: f64) -> f64 {
        0.5 * self.inductance() * current * current
    }

    pub fn magnetic_field(&self, current: f64) -> f64 {
        MU_0 * self.core_mu_r * (self.turns / self.length) * current
    }

    pub fn mutual_inductance(&self, other_turns: f64) -> f64 {
        MU_0 * self.core_mu_r * self.turns * other_turns * self.area / self.length
    }
}

/// Capacitors in series: 1/C = Σ 1/Cᵢ
pub fn capacitors_series(capacitances: &[f64]) -> f64 {
    1.0 / capacitances.iter().map(|c| 1.0 / c).sum::<f64>()
}

/// Capacitors in parallel: C = Σ Cᵢ
pub fn capacitors_parallel(capacitances: &[f64]) -> f64 {
    capacitances.iter().sum()
}

/// Inductors in series
pub fn inductors_series(inductances: &[f64]) -> f64 {
    inductances.iter().sum()
}

/// Inductors in parallel
pub fn inductors_parallel(inductances: &[f64]) -> f64 {
    1.0 / inductances.iter().map(|l| 1.0 / l).sum::<f64>()
}

/// Energy in capacitor: U = ½CV²
pub fn capacitor_energy(capacitance: f64, voltage: f64) -> f64 {
    0.5 * capacitance * voltage * voltage
}

/// Energy in inductor: U = ½LI²
pub fn inductor_energy(inductance: f64, current: f64) -> f64 {
    0.5 * inductance * current * current
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parallel_plate_capacitance() {
        let cap = ParallelPlateCapacitor::new(1.0, 0.01);
        let expected = EPSILON_0 / 0.01;
        assert!((cap.capacitance() - expected).abs() / expected < 1e-10);
    }

    #[test]
    fn test_capacitor_with_dielectric() {
        let cap = ParallelPlateCapacitor::with_dielectric(1.0, 0.01, 4.0);
        let expected = EPSILON_0 * 4.0 / 0.01;
        assert!((cap.capacitance() - expected).abs() / expected < 1e-10);
    }

    #[test]
    fn test_capacitor_energy() {
        let cap = ParallelPlateCapacitor::new(1.0, 0.01);
        let u = cap.energy(10.0);
        assert!((u - 0.5 * cap.capacitance() * 100.0).abs() < 1e-15);
    }

    #[test]
    fn test_spherical_capacitor() {
        let cap = SphericalCapacitor::new(0.05, 0.1);
        let expected = 4.0 * std::f64::consts::PI * EPSILON_0 * 0.05 * 0.1 / 0.05;
        assert!((cap.capacitance() - expected).abs() / expected < 1e-10);
    }

    #[test]
    fn test_solenoid_inductance() {
        let sol = SolenoidInductor::new(100.0, 1.0, 0.01);
        let expected = MU_0 * 10000.0 * 0.01;
        assert!((sol.inductance() - expected).abs() / expected < 1e-10);
    }

    #[test]
    fn test_solenoid_with_core() {
        let sol = SolenoidInductor::with_core(100.0, 1.0, 0.01, 100.0);
        let expected = MU_0 * 100.0 * 10000.0 * 0.01;
        assert!((sol.inductance() - expected).abs() / expected < 1e-10);
    }

    #[test]
    fn test_capacitors_series() {
        assert!((capacitors_series(&[1e-6, 1e-6]) - 0.5e-6).abs() < 1e-15);
    }

    #[test]
    fn test_capacitors_parallel() {
        assert!((capacitors_parallel(&[1e-6, 2e-6]) - 3e-6).abs() < 1e-15);
    }

    #[test]
    fn test_inductors_series() {
        assert!((inductors_series(&[1e-3, 2e-3]) - 3e-3).abs() < 1e-15);
    }

    #[test]
    fn test_inductors_parallel() {
        assert!((inductors_parallel(&[1e-3, 1e-3]) - 0.5e-3).abs() < 1e-15);
    }
}
