use crate::constants::MU_0;
use crate::vector::Vec3;

/// A current-carrying wire element for Biot-Savart computation.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct WireElement {
    pub start: Vec3,
    pub end: Vec3,
    pub current: f64,
}

impl WireElement {
    pub fn new(start: Vec3, end: Vec3, current: f64) -> Self {
        Self { start, end, current }
    }

    pub fn direction(&self) -> Vec3 {
        let d = self.end - self.start;
        if d.norm() < f64::EPSILON {
            Vec3::zeros()
        } else {
            d.normalize()
        }
    }

    pub fn length(&self) -> f64 {
        (self.end - self.start).norm()
    }

    pub fn midpoint(&self) -> Vec3 {
        (self.start + self.end) / 2.0
    }
}

/// Biot-Savart law: compute magnetic field at `point` from a wire element.
pub fn biot_savart_segment(wire: &WireElement, point: Vec3) -> Vec3 {
    let dl = wire.end - wire.start;
    let n_segments = 100;
    let mut total_field = Vec3::zeros();

    for i in 0..n_segments {
        let t0 = i as f64 / n_segments as f64;
        let t1 = (i + 1) as f64 / n_segments as f64;
        let seg_start = wire.start + t0 * dl;
        let seg_end = wire.start + t1 * dl;
        let seg_dl = seg_end - seg_start;
        let seg_mid = (seg_start + seg_end) / 2.0;

        let r_vec = point - seg_mid;
        let r = r_vec.norm();
        if r < f64::EPSILON {
            continue;
        }
        let r_hat = r_vec / r;
        let cross = seg_dl.cross(&r_hat);
        total_field += (MU_0 / (4.0 * std::f64::consts::PI)) * wire.current * cross / (r * r);
    }

    total_field
}

/// Magnetic field from an infinite straight wire at perpendicular distance d.
/// B = μ₀ I / (2π d)
pub fn infinite_wire_field(current: f64, distance: f64) -> f64 {
    MU_0 * current / (2.0 * std::f64::consts::PI * distance)
}

/// Magnetic field at center of a circular current loop.
/// B = μ₀ I / (2R)
pub fn loop_center_field(current: f64, radius: f64) -> f64 {
    MU_0 * current / (2.0 * radius)
}

/// Ampere's law: ∮ B·dl = μ₀ I_enc
pub fn amperes_law_circular(current: f64, _radius: f64) -> f64 {
    MU_0 * current
}

/// Magnetic field from a magnetic dipole on axis.
pub fn magnetic_dipole_axial(moment: f64, z: f64) -> f64 {
    if z.abs() < f64::EPSILON {
        return f64::INFINITY;
    }
    (MU_0 / (4.0 * std::f64::consts::PI)) * 2.0 * moment / (z * z * z)
}

/// Solenoid field: B = μ₀ n I
pub fn solenoid_field(turns_per_meter: f64, current: f64) -> f64 {
    MU_0 * turns_per_meter * current
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::vector;

    #[test]
    fn test_infinite_wire_field() {
        let b = infinite_wire_field(1.0, 1.0);
        let expected = MU_0 / (2.0 * std::f64::consts::PI);
        assert!((b - expected).abs() / expected < 1e-10);
    }

    #[test]
    fn test_infinite_wire_inverse_distance() {
        let b1 = infinite_wire_field(1.0, 1.0);
        let b2 = infinite_wire_field(1.0, 2.0);
        assert!((b1 / b2 - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_loop_center_field() {
        let b = loop_center_field(1.0, 0.1);
        let expected = MU_0 * 1.0 / (2.0 * 0.1);
        assert!((b - expected).abs() / expected < 1e-10);
    }

    #[test]
    fn test_amperes_law() {
        let circulation = amperes_law_circular(5.0, 0.5);
        let expected = MU_0 * 5.0;
        assert!((circulation - expected).abs() / expected < 1e-10);
    }

    #[test]
    fn test_solenoid_field() {
        let b = solenoid_field(1000.0, 1.0);
        let expected = MU_0 * 1000.0;
        assert!((b - expected).abs() / expected < 1e-10);
    }

    #[test]
    fn test_biot_savart_straight_wire() {
        let wire = WireElement::new(
            vector![0.0, 0.0, -1.0],
            vector![0.0, 0.0, 1.0],
            1.0,
        );
        let point = vector![1.0, 0.0, 0.0];
        let b = biot_savart_segment(&wire, point);
        assert!(b.norm() > 0.0);
    }

    #[test]
    fn test_magnetic_dipole_cube_law() {
        let b1 = magnetic_dipole_axial(1.0, 1.0);
        let b2 = magnetic_dipole_axial(1.0, 2.0);
        assert!((b1 / b2 - 8.0).abs() < 1e-6);
    }
}
