use crate::constants::{C, EPSILON_0, MU_0};
use crate::vector::Vec3;

/// Electromagnetic wave parameters.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct EMWave {
    pub frequency: f64,
    pub wavelength: f64,
    pub amplitude_e: f64,
    pub amplitude_b: f64,
}

impl EMWave {
    pub fn from_frequency(f: f64, amplitude_e: f64) -> Self {
        let wavelength = C / f;
        let amplitude_b = amplitude_e / C;
        Self { frequency: f, wavelength, amplitude_e, amplitude_b }
    }

    pub fn from_wavelength(lambda: f64, amplitude_e: f64) -> Self {
        let frequency = C / lambda;
        let amplitude_b = amplitude_e / C;
        Self { frequency, wavelength: lambda, amplitude_e, amplitude_b }
    }

    pub fn angular_frequency(&self) -> f64 {
        2.0 * std::f64::consts::PI * self.frequency
    }

    pub fn wave_number(&self) -> f64 {
        2.0 * std::f64::consts::PI / self.wavelength
    }

    pub fn period(&self) -> f64 {
        1.0 / self.frequency
    }

    pub fn verify_eb_ratio(&self) -> bool {
        if self.amplitude_b.abs() < f64::EPSILON {
            return false;
        }
        let ratio = self.amplitude_e / self.amplitude_b;
        (ratio - C).abs() / C < 1e-6
    }

    pub fn e_field(&self, x: f64, t: f64, polarization: &Polarization) -> Vec3 {
        let phase = self.wave_number() * x - self.angular_frequency() * t;
        match polarization {
            Polarization::LinearY => Vec3::new(0.0, self.amplitude_e * phase.cos(), 0.0),
            Polarization::LinearZ => Vec3::new(0.0, 0.0, self.amplitude_e * phase.cos()),
            Polarization::Circular(direction) => {
                let dir = if *direction == Handedness::Right { 1.0 } else { -1.0 };
                Vec3::new(0.0, self.amplitude_e * phase.cos(), dir * self.amplitude_e * phase.sin())
            }
        }
    }

    pub fn b_field(&self, x: f64, t: f64, polarization: &Polarization) -> Vec3 {
        let phase = self.wave_number() * x - self.angular_frequency() * t;
        match polarization {
            Polarization::LinearY => Vec3::new(0.0, 0.0, self.amplitude_b * phase.cos()),
            Polarization::LinearZ => Vec3::new(0.0, -self.amplitude_b * phase.cos(), 0.0),
            Polarization::Circular(direction) => {
                let dir = if *direction == Handedness::Right { 1.0 } else { -1.0 };
                Vec3::new(0.0, -dir * self.amplitude_b * phase.sin(), self.amplitude_b * phase.cos())
            }
        }
    }
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub enum Polarization {
    LinearY,
    LinearZ,
    Circular(Handedness),
}

#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum Handedness {
    Right,
    Left,
}

/// Poynting vector: S = (1/μ₀) E × B
pub fn poynting_vector(e: &Vec3, b: &Vec3) -> Vec3 {
    e.cross(b) / MU_0
}

/// Time-averaged intensity: I = E₀² / (2μ₀c)
pub fn average_intensity(amplitude_e: f64) -> f64 {
    amplitude_e * amplitude_e / (2.0 * MU_0 * C)
}

/// Energy density of electric field: u_E = ½ε₀E²
pub fn electric_energy_density(e_magnitude: f64) -> f64 {
    0.5 * EPSILON_0 * e_magnitude * e_magnitude
}

/// Energy density of magnetic field: u_B = B²/(2μ₀)
pub fn magnetic_energy_density(b_magnitude: f64) -> f64 {
    b_magnitude * b_magnitude / (2.0 * MU_0)
}

/// Total EM energy density
pub fn total_energy_density(e: f64, b: f64) -> f64 {
    electric_energy_density(e) + magnetic_energy_density(b)
}

/// Radiation pressure on absorbing surface
pub fn radiation_pressure_absorb(intensity: f64) -> f64 {
    intensity / C
}

/// Radiation pressure on reflecting surface
pub fn radiation_pressure_reflect(intensity: f64) -> f64 {
    2.0 * intensity / C
}

/// Verify u_E = u_B in an EM wave
pub fn verify_equal_energy_densities(e: f64, b: f64) -> bool {
    let u_e = electric_energy_density(e);
    let u_b = magnetic_energy_density(b);
    if u_e.abs() < f64::EPSILON {
        return u_b.abs() < 1e-10;
    }
    (u_e - u_b).abs() / u_e < 1e-6
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::vector;

    #[test]
    fn test_wave_frequency_wavelength() {
        let wave = EMWave::from_frequency(1e9, 1.0);
        assert!((wave.wavelength - 0.299_792_458).abs() < 1e-3);
    }

    #[test]
    fn test_wave_from_wavelength() {
        let wave = EMWave::from_wavelength(0.5, 1.0);
        let expected_freq = C / 0.5;
        assert!((wave.frequency - expected_freq).abs() / expected_freq < 1e-6);
    }

    #[test]
    fn test_eb_ratio() {
        let wave = EMWave::from_frequency(1e9, 100.0);
        assert!(wave.verify_eb_ratio());
    }

    #[test]
    fn test_poynting_vector() {
        let e = vector![0.0, 1.0, 0.0];
        let b = vector![0.0, 0.0, 1.0 / C];
        let s = poynting_vector(&e, &b);
        assert!(s.x > 0.0);
        assert!(s.y.abs() < 1e-15);
        assert!(s.z.abs() < 1e-15);
    }

    #[test]
    fn test_poynting_magnitude() {
        let e = vector![0.0, 100.0, 0.0];
        let b = vector![0.0, 0.0, 100.0 / C];
        let s = poynting_vector(&e, &b);
        let expected = 100.0 * (100.0 / C) / MU_0;
        assert!((s.norm() - expected).abs() / expected < 1e-6);
    }

    #[test]
    fn test_average_intensity() {
        let i = average_intensity(100.0);
        let expected = 100.0 * 100.0 / (2.0 * MU_0 * C);
        assert!((i - expected).abs() / expected < 1e-6);
    }

    #[test]
    fn test_energy_density_equality() {
        let e = 100.0;
        let b = e / C;
        assert!(verify_equal_energy_densities(e, b));
    }

    #[test]
    fn test_electric_energy_density() {
        let u = electric_energy_density(1.0);
        let expected = 0.5 * EPSILON_0;
        assert!((u - expected).abs() < 1e-15);
    }

    #[test]
    fn test_radiation_pressure_absorb() {
        let i = average_intensity(100.0);
        let p = radiation_pressure_absorb(i);
        assert!(p > 0.0);
        assert!((p - i / C).abs() < 1e-15);
    }

    #[test]
    fn test_radiation_pressure_reflect() {
        let p = radiation_pressure_reflect(1.0);
        assert!((p - 2.0 / C).abs() < 1e-15);
    }

    #[test]
    fn test_plane_wave_propagation() {
        let wave = EMWave::from_frequency(1e9, 1.0);
        let pol = Polarization::LinearY;
        let e0 = wave.e_field(0.0, 0.0, &pol);
        assert!((e0.y - 1.0).abs() < 1e-10);
        let e1 = wave.e_field(wave.wavelength, 0.0, &pol);
        assert!((e1.y - 1.0).abs() < 1e-6);
    }
}
