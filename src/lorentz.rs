use crate::constants::{E_CHARGE, ELECTRON_MASS, PROTON_MASS};
use crate::vector::Vec3;

/// Lorentz force: F = q(E + v × B)
pub fn lorentz_force(charge: f64, e_field: &Vec3, b_field: &Vec3, velocity: &Vec3) -> Vec3 {
    charge * (e_field + velocity.cross(b_field))
}

/// Electric force: F = qE
pub fn electric_force(charge: f64, e_field: &Vec3) -> Vec3 {
    charge * e_field
}

/// Magnetic force: F = qv × B
pub fn magnetic_force(charge: f64, velocity: &Vec3, b_field: &Vec3) -> Vec3 {
    charge * velocity.cross(b_field)
}

/// Cyclotron radius: r = mv/(|q|B)
pub fn cyclotron_radius(mass: f64, speed: f64, charge: f64, b_magnitude: f64) -> f64 {
    if charge.abs() < f64::EPSILON || b_magnitude.abs() < f64::EPSILON {
        return f64::INFINITY;
    }
    mass * speed / (charge.abs() * b_magnitude)
}

/// Cyclotron frequency: ω = |q|B/m
pub fn cyclotron_frequency(charge: f64, b_magnitude: f64, mass: f64) -> f64 {
    if mass.abs() < f64::EPSILON {
        return f64::INFINITY;
    }
    charge.abs() * b_magnitude / mass
}

/// Cyclotron period: T = 2π/ω
pub fn cyclotron_period(charge: f64, b_magnitude: f64, mass: f64) -> f64 {
    2.0 * std::f64::consts::PI / cyclotron_frequency(charge, b_magnitude, mass)
}

/// Charged particle state.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Particle {
    pub charge: f64,
    pub mass: f64,
    pub position: Vec3,
    pub velocity: Vec3,
}

impl Particle {
    pub fn new(charge: f64, mass: f64, position: Vec3, velocity: Vec3) -> Self {
        Self { charge, mass, position, velocity }
    }

    pub fn electron(position: Vec3, velocity: Vec3) -> Self {
        Self::new(-E_CHARGE, ELECTRON_MASS, position, velocity)
    }

    pub fn proton(position: Vec3, velocity: Vec3) -> Self {
        Self::new(E_CHARGE, PROTON_MASS, position, velocity)
    }

    /// Boris algorithm step for magnetic field integration (preserves energy).
    pub fn step_boris(&self, e_field: &Vec3, b_field: &Vec3, dt: f64) -> Particle {
        let qm = self.charge / self.mass;
        let v_minus = self.velocity + 0.5 * qm * e_field * dt;
        let t_vec = 0.5 * qm * b_field * dt;
        let t_sq = t_vec.dot(&t_vec);
        let s_vec = 2.0 * t_vec / (1.0 + t_sq);
        let v_prime = v_minus + v_minus.cross(&t_vec);
        let v_plus = v_minus + v_prime.cross(&s_vec);
        let new_velocity = v_plus + 0.5 * qm * e_field * dt;
        let new_position = self.position + new_velocity * dt;
        Particle { charge: self.charge, mass: self.mass, position: new_position, velocity: new_velocity }
    }

    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * self.velocity.norm_squared()
    }

    pub fn speed(&self) -> f64 {
        self.velocity.norm()
    }
}

/// Helical motion parameters.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct HelicalMotion {
    pub radius: f64,
    pub frequency: f64,
    pub pitch: f64,
}

impl HelicalMotion {
    pub fn compute(particle: &Particle, b_magnitude: f64, v_perp: f64, v_parallel: f64) -> Self {
        let r = cyclotron_radius(particle.mass, v_perp, particle.charge, b_magnitude);
        let f = cyclotron_frequency(particle.charge, b_magnitude, particle.mass);
        let pitch = v_parallel / f;
        Self { radius: r, frequency: f, pitch }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::vector;

    #[test]
    fn test_lorentz_force_electric() {
        let f = lorentz_force(1.0, &vector![1.0, 0.0, 0.0], &Vec3::zeros(), &Vec3::zeros());
        assert!((f.x - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_lorentz_force_magnetic() {
        let v = vector![1.0, 0.0, 0.0];
        let b = vector![0.0, 0.0, 1.0];
        let f = lorentz_force(1.0, &Vec3::zeros(), &b, &v);
        assert!(f.x.abs() < 1e-10);
        assert!((f.y - (-1.0)).abs() < 1e-10);
        assert!(f.z.abs() < 1e-10);
    }

    #[test]
    fn test_magnetic_force_perpendicular() {
        let v = vector![1.0, 0.0, 0.0];
        let b = vector![0.0, 0.0, 1.0];
        let f = magnetic_force(1.0, &v, &b);
        assert!(f.dot(&v).abs() < 1e-10);
    }

    #[test]
    fn test_cyclotron_radius() {
        assert!((cyclotron_radius(1.0, 1.0, 1.0, 1.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_cyclotron_frequency() {
        assert!((cyclotron_frequency(1.0, 1.0, 1.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_cyclotron_period() {
        assert!((cyclotron_period(1.0, 1.0, 1.0) - 2.0 * std::f64::consts::PI).abs() < 1e-10);
    }

    #[test]
    fn test_boris_algorithm_energy_conservation() {
        let particle = Particle::new(1.0, 1.0, vector![0.0, 0.0, 0.0], vector![1.0, 0.0, 0.0]);
        let b = vector![0.0, 0.0, 1.0];
        let initial_ke = particle.kinetic_energy();
        let mut p = particle.clone();
        for _ in 0..100 {
            p = p.step_boris(&Vec3::zeros(), &b, 0.01);
        }
        assert!((p.kinetic_energy() - initial_ke).abs() / initial_ke < 0.01);
    }

    #[test]
    fn test_helical_motion() {
        let particle = Particle::new(1.0, 1.0, Vec3::zeros(), vector![1.0, 0.0, 1.0]);
        let helix = HelicalMotion::compute(&particle, 1.0, 1.0, 1.0);
        assert!(helix.radius > 0.0);
        assert!(helix.frequency > 0.0);
        assert!(helix.pitch > 0.0);
    }

    #[test]
    fn test_particle_electron_creation() {
        let e = Particle::electron(Vec3::zeros(), vector![1e6, 0.0, 0.0]);
        assert!((e.charge - (-E_CHARGE)).abs() < 1e-30);
        assert!((e.mass - ELECTRON_MASS).abs() < 1e-40);
    }
}
