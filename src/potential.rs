use crate::constants::{EPSILON_0, K_COULOMB};
use crate::vector::{distance, Vec3};

/// Electric scalar potential from a point charge: V = kQ/r
pub fn scalar_potential_point_charge(charge: f64, r: f64) -> f64 {
    if r.abs() < f64::EPSILON {
        return f64::INFINITY.signum() * charge;
    }
    K_COULOMB * charge / r
}

/// Potential from multiple charges (superposition).
pub fn scalar_potential_superposition(charges: &[(Vec3, f64)], point: Vec3) -> f64 {
    charges.iter().map(|(pos, q)| scalar_potential_point_charge(*q, distance(*pos, point))).sum()
}

/// Potential difference: ΔV = V(b) - V(a)
pub fn potential_difference(v_b: f64, v_a: f64) -> f64 {
    v_b - v_a
}

/// E = -dV/dx
pub fn e_field_from_potential_1d(dv_dx: f64) -> f64 {
    -dv_dx
}

/// Poisson's equation residual: ∇²V + ρ/ε₀ (should be zero)
pub fn poisson_equation(laplacian_v: f64, rho: f64) -> f64 {
    laplacian_v + rho / EPSILON_0
}

/// Laplace's equation residual (should be zero)
pub fn laplace_equation(laplacian_v: f64) -> f64 {
    laplacian_v
}

/// Potential of a uniformly charged sphere at distance r.
pub fn potential_charged_sphere(q: f64, r: f64, radius: f64) -> f64 {
    if r >= radius {
        K_COULOMB * q / r
    } else {
        K_COULOMB * q * (3.0 * radius * radius - r * r) / (2.0 * radius * radius * radius)
    }
}

/// Electric potential energy: U = kq₁q₂/r
pub fn potential_energy_pair(q1: f64, q2: f64, r: f64) -> f64 {
    K_COULOMB * q1 * q2 / r
}

/// Work from potential: W = qΔV
pub fn work_from_potential(charge: f64, delta_v: f64) -> f64 {
    charge * delta_v
}

/// Numerical gradient of scalar potential (3D central differences).
pub fn potential_gradient(v_func: &dyn Fn(Vec3) -> f64, point: Vec3, h: f64) -> Vec3 {
    use nalgebra::vector;
    let dx = vector![h, 0.0, 0.0];
    let dy = vector![0.0, h, 0.0];
    let dz = vector![0.0, 0.0, h];
    Vec3::new(
        (v_func(point + dx) - v_func(point - dx)) / (2.0 * h),
        (v_func(point + dy) - v_func(point - dy)) / (2.0 * h),
        (v_func(point + dz) - v_func(point - dz)) / (2.0 * h),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::vector;

    #[test]
    fn test_point_charge_potential() {
        let v = scalar_potential_point_charge(1e-6, 1.0);
        let expected = K_COULOMB * 1e-6;
        assert!((v - expected).abs() / expected < 1e-10);
    }

    #[test]
    fn test_potential_superposition() {
        let charges = vec![
            (vector![1.0, 0.0, 0.0], 1e-6),
            (vector![-1.0, 0.0, 0.0], -1e-6),
        ];
        let v = scalar_potential_superposition(&charges, vector![0.0, 0.0, 0.0]);
        assert!(v.abs() < 1e-10);
    }

    #[test]
    fn test_potential_difference() {
        assert!((potential_difference(10.0, 5.0) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_e_field_from_potential() {
        assert!((e_field_from_potential_1d(-100.0) - 100.0).abs() < 1e-10);
    }

    #[test]
    fn test_poisson_equation() {
        let rho = 1e-6;
        let residual = poisson_equation(-rho / EPSILON_0, rho);
        assert!(residual.abs() < 1e-10);
    }

    #[test]
    fn test_laplace_equation() {
        assert!(laplace_equation(0.0).abs() < 1e-10);
    }

    #[test]
    fn test_charged_sphere_outside() {
        let v = potential_charged_sphere(1e-6, 2.0, 1.0);
        let expected = K_COULOMB * 1e-6 / 2.0;
        assert!((v - expected).abs() / expected < 1e-10);
    }

    #[test]
    fn test_charged_sphere_inside() {
        let v = potential_charged_sphere(1e-6, 0.0, 1.0);
        let v_surface = potential_charged_sphere(1e-6, 1.0, 1.0);
        assert!(v > v_surface);
    }

    #[test]
    fn test_potential_energy() {
        let u = potential_energy_pair(1e-6, 1e-6, 1.0);
        let expected = K_COULOMB * 1e-12;
        assert!((u - expected).abs() / expected < 1e-10);
    }

    #[test]
    fn test_work_from_potential() {
        assert!((work_from_potential(1e-6, 100.0) - 1e-4).abs() < 1e-15);
    }

    #[test]
    fn test_potential_gradient() {
        let v_func = |p: Vec3| {
            let r = p.norm();
            if r < 1e-10 { 0.0 } else { K_COULOMB * 1e-6 / r }
        };
        let grad = potential_gradient(&v_func, vector![1.0, 0.0, 0.0], 1e-6);
        let e = -grad;
        assert!(e.x > 0.0);
    }
}
