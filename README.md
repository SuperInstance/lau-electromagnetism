# lau-electromagnetism

A classical electromagnetism library in Rust. Maxwell's equations, Coulomb's law, Biot-Savart, Lorentz force, EM waves, potentials, capacitance, materials, and a novel **agent-field interaction** model — all built on `nalgebra` 3-vectors with full serde support.

**90 tests** · MIT licensed · `nalgebra` + `serde` only.

---

## Table of Contents

1. [What This Does](#what-this-does)
2. [Key Idea](#key-idea)
3. [Install](#install)
4. [Quick Start](#quick-start)
5. [API Reference](#api-reference)
6. [How It Works](#how-it-works)
7. [The Math](#the-math)
8. [License](#license)

---

## What This Does

The crate is organised into 11 modules:

| Module | What it covers |
|---|---|
| `vector` | 3D vector helpers (`Vec3`, distance, unit vector, cross, dot) |
| `constants` | ε₀, μ₀, c, k_Coulomb, e, m_e, m_p, k_B, Z₀ |
| `electrostatics` | Coulomb force, electric field (point + superposition), Gauss's law, line/plane charges |
| `magnetostatics` | Biot-Savart, infinite wire, current loop, Ampère's law, solenoid, magnetic dipole |
| `maxwell` | Maxwell's equations (differential + integral), displacement current, induced EMF, magnetic flux |
| `lorentz` | Lorentz force, cyclotron motion, Boris integrator, helical motion |
| `potential` | Scalar potential, superposition, Poisson/Laplace, charged sphere, numerical gradient |
| `capacitance` | Parallel plate, spherical, cylindrical capacitors; solenoid inductor; series/parallel combinations |
| `materials` | Permittivity, permeability, skin depth, refractive index, boundary conditions; built-in materials (vacuum, glass, copper, …) |
| `waves` | EM wave parameters, Poynting vector, energy density, radiation pressure, polarisation |
| `agent_field` | Agent-as-charge model: simulate N charged agents interacting via Coulomb fields |

---

## Key Idea

This library treats classical EM as a **compute graph**: you compose point charges, wire elements, capacitor geometries, and material properties, then call functions that return exact analytical results (or numerical approximations where analytical solutions don't exist, e.g. Biot-Savart integration).

The `agent_field` module extends this into an agent-based modelling paradigm: each agent carries an "influence charge" and a sensitivity, and agents interact via the same Coulomb fields — repulsion, attraction, and superposition all work. This lets you simulate crowd dynamics, social influence spreading, or any system that maps onto inverse-square-law interactions.

---

## Install

```toml
[dependencies]
lau-electromagnetism = "0.1"
```

### Dependencies

| Crate | Why |
|---|---|
| `nalgebra` 0.33 (with `serde-serialize`) | 3D vector math |
| `serde` 1 (with `derive`) | Serialisation |

---

## Quick Start

### Electrostatics

```rust
use lau_electromagnetism::*;
use nalgebra::vector;

let q1 = PointCharge::new(vector![0.0, 0.0, 0.0], 1e-6);
let q2 = PointCharge::new(vector![1.0, 0.0, 0.0], -1e-6);

// Coulomb force on q2 due to q1
let force = coulomb_force(&q1, &q2);
assert!(force.x < 0.0); // opposite charges attract

// Electric field at a point from multiple charges
let field = electric_field_superposition(&[q1, q2], vector![0.5, 0.0, 0.0]);

// Gauss's law: flux through a sphere
let flux = electric_flux_sphere(&charges, 1.0);
```

### Magnetostatics

```rust
let wire = WireElement::new(vector![0.0, 0.0, -1.0], vector![0.0, 0.0, 1.0], 1.0);
let b_field = biot_savart_segment(&wire, vector![1.0, 0.0, 0.0]);

// Infinite wire: B = μ₀I / (2πd)
let b = infinite_wire_field(1.0, 0.1);

// Solenoid: B = μ₀nI
let b = solenoid_field(1000.0, 1.0);
```

### Lorentz Force & Boris Integrator

```rust
let electron = Particle::electron(
    vector![0.0, 0.0, 0.0],
    vector![1e6, 0.0, 0.0],
);
let b = vector![0.0, 0.0, 1.0];

// Boris pusher (energy-conserving)
let mut p = electron;
for _ in 0..100 {
    p = p.step_boris(&Vec3::zeros(), &b, 1e-9);
}
// Kinetic energy preserved to < 1%
```

### Agent Field Simulation

```rust
let mut agents = vec![
    Agent::new("a1", vector![-1.0, 0.0, 0.0], 1.0, 1.0),
    Agent::new("a2", vector![ 1.0, 0.0, 0.0], 1.0, 1.0),
];

// Like charges repel — agents drift apart
simulate_agent_interactions(&mut agents, 0.001, 100);
```

---

## API Reference

### `vector` module

| Function | Signature | Returns |
|---|---|---|
| `unit_vector(from, to)` | `Vec3, Vec3 → Vec3` | Normalised direction, or zero if coincident |
| `distance(a, b)` | `Vec3, Vec3 → f64` | Euclidean distance |
| `cross(a, b)` | `&Vec3, &Vec3 → Vec3` | Cross product |
| `dot(a, b)` | `&Vec3, &Vec3 → f64` | Dot product |

### `constants` module

| Constant | Value | Units |
|---|---|---|
| `EPSILON_0` | 8.854187817 × 10⁻¹² | F/m |
| `MU_0` | 4π × 10⁻⁷ | H/m |
| `C` | 299,792,458 | m/s |
| `K_COULOMB` | 1/(4πε₀) ≈ 8.988 × 10⁹ | N·m²/C² |
| `E_CHARGE` | 1.602176634 × 10⁻¹⁹ | C |
| `ELECTRON_MASS` | 9.1093837015 × 10⁻³¹ | kg |
| `PROTON_MASS` | 1.67262192369 × 10⁻²⁷ | kg |
| `BOLTZMANN` | 1.380649 × 10⁻²³ | J/K |
| `Z_0` | μ₀c ≈ 376.73 | Ω |

`verify_speed_of_light()` checks that `1/√(ε₀μ₀) = c` within tolerance.

### `electrostatics` module

| Function | Formula |
|---|---|
| `coulomb_force(q1, q2)` | F = k q₁q₂ / r² r̂ |
| `electric_field_point(charge, point)` | E = k q / r² r̂ |
| `electric_field_superposition(charges, point)` | Σ Eᵢ |
| `electric_flux_sphere(charges, radius)` | Φ = Q_enc / ε₀ |
| `electric_field_line_charge(λ, d)` | E = λ / (2πε₀d) |
| `electric_field_plane_charge(σ)` | E = σ / (2ε₀) |

### `magnetostatics` module

| Function | Formula |
|---|---|
| `biot_savart_segment(wire, point)` | Numerical integration (100 segments) |
| `infinite_wire_field(I, d)` | B = μ₀I / (2πd) |
| `loop_center_field(I, R)` | B = μ₀I / (2R) |
| `solenoid_field(n, I)` | B = μ₀nI |
| `magnetic_dipole_axial(m, z)` | B = (μ₀/4π) 2m / z³ |

### `maxwell` module

`MaxwellFields` struct holds E, B, J, ρ, ∂E/∂t, ∂B/∂t and provides residual checks for all four Maxwell's equations:

| Method | Equation |
|---|---|
| `gauss_law_divergence(div_e)` | ∇·E − ρ/ε₀ = 0 |
| `gauss_magnetism(div_b)` | ∇·B = 0 |
| `faradays_law(curl_e)` | ∇×E + ∂B/∂t = 0 |
| `ampere_maxwell_law(curl_b)` | ∇×B − μ₀J − μ₀ε₀∂E/∂t = 0 |

Integral forms: `gauss_law_integral`, `faradays_law_integral`, `amperes_law_integral`.

### `lorentz` module

| Function/Type | Description |
|---|---|
| `lorentz_force(q, E, B, v)` | F = q(E + v × B) |
| `cyclotron_radius(m, v, q, B)` | r = mv / (\|q\|B) |
| `cyclotron_frequency(q, B, m)` | ω = \|q\|B / m |
| `Particle::step_boris(E, B, dt)` | Boris pusher (symplectic, energy-conserving) |
| `Particle::electron(pos, vel)` | Convenience constructor |
| `HelicalMotion::compute(...)` | Radius, frequency, pitch of helical orbit |

### `potential` module

| Function | Formula |
|---|---|
| `scalar_potential_point_charge(q, r)` | V = kq/r |
| `scalar_potential_superposition(charges, point)` | Σ Vᵢ |
| `potential_charged_sphere(Q, r, R)` | Inside/outside analytic formula |
| `potential_energy_pair(q₁, q₂, r)` | U = kq₁q₂/r |
| `potential_gradient(v_func, point, h)` | Numerical ∂V/∂x via central differences |
| `poisson_equation(∇²V, ρ)` | Residual check: ∇²V + ρ/ε₀ = 0 |

### `capacitance` module

| Type | Formula |
|---|---|
| `ParallelPlateCapacitor` | C = ε₀κA/d |
| `SphericalCapacitor` | C = 4πε₀κ ab/(b−a) |
| `CylindricalCapacitor` | C = 2πε₀κL / ln(b/a) |
| `SolenoidInductor` | L = μ₀μᵣN²A/ℓ |

Plus `capacitors_series`, `capacitors_parallel`, `inductors_series`, `inductors_parallel`.

### `materials` module

Pre-built materials: `vacuum()`, `air()`, `glass()`, `water()`, `copper()`, `iron()`, `teflon()`.

Each `Material` provides: `permittivity()`, `permeability()`, `wave_speed()`, `refractive_index()`, `intrinsic_impedance()`, `skin_depth(f)`, `wavelength_in_material(λ₀)`.

### `waves` module

| Function | Formula |
|---|---|
| `EMWave::from_frequency(f, E₀)` | λ = c/f, B₀ = E₀/c |
| `poynting_vector(E, B)` | S = (1/μ₀) E × B |
| `average_intensity(E₀)` | I = E₀² / (2μ₀c) |
| `electric_energy_density(E)` | u_E = ½ε₀E² |
| `magnetic_energy_density(B)` | u_B = B²/(2μ₀) |
| `radiation_pressure_absorb(I)` | P = I/c |
| `radiation_pressure_reflect(I)` | P = 2I/c |

Supports `LinearY`, `LinearZ`, and `Circular(Handedness)` polarisation.

### `agent_field` module

| Type/Function | Description |
|---|---|
| `Agent` | id, position, influence_charge, sensitivity, velocity |
| `FieldMap` | Grid of agents with field/potential sampling |
| `simulate_agent_interactions(agents, dt, steps)` | N-body Coulomb simulation |
| `interaction_energy(a, b)` | U = kq₁q₂/r |
| `total_interaction_energy(agents)` | Sum of all pairwise energies |

---

## How It Works

### Design Philosophy

The library separates **pure functions** (given inputs, compute outputs) from **stateful types** (particles, agents, capacitors). This makes everything testable: every formula is checked against known analytical results.

### Biot-Savart Numerical Integration

For arbitrary wire segments, the library subdivides each `WireElement` into 100 straight segments and sums:

```
B_total = Σ (μ₀ / 4π) · I · (dl × r̂) / r²
```

This converges well for smooth geometries. For the special cases (infinite wire, loop, solenoid), closed-form solutions are provided instead.

### Boris Pusher

The `Particle::step_boris` method implements the **Boris algorithm**, the gold-standard integrator for charged particle motion in magnetic fields. It splits the velocity update into:

1. Half-kick from E field
2. Rotation from B field (exactly preserves |v|)
3. Half-kick from E field

This is **symplectic**: kinetic energy is conserved to machine precision over long runs, unlike naive Euler which spirals outward.

### Agent-Field Analogy

The `agent_field` module maps each agent to a point charge. Agent–agent forces use the same `coulomb_force` function. This creates an N-body problem solved by direct summation (O(n²) per step). The `simulate_agent_interactions` function is a simple leapfrog integrator.

---

## The Math

### Maxwell's Equations (differential form)

| Equation | Form |
|---|---|
| Gauss (E) | ∇ · **E** = ρ / ε₀ |
| Gauss (B) | ∇ · **B** = 0 |
| Faraday | ∇ × **E** = −∂**B**/∂t |
| Ampère-Maxwell | ∇ × **B** = μ₀**J** + μ₀ε₀ ∂**E**/∂t |

### Coulomb's Law

```
F₁₂ = k · q₁q₂ / r² · r̂₁₂
```

where k = 1/(4πε₀) ≈ 8.988 × 10⁹ N·m²/C².

### Lorentz Force

```
F = q(E + v × B)
```

The magnetic force `qv × B` is always perpendicular to **v**, so it does no work — it changes direction but not speed.

### Cyclotron Motion

A charged particle in a uniform B field traces a circle (or helix if there's a velocity component along B):

```
r = mv / (|q|B)        (cyclotron radius)
ω = |q|B / m           (cyclotron frequency)
T = 2πm / (|q|B)       (cyclotron period)
```

### Energy in Fields

```
u_E = ½ε₀E²            (electric energy density)
u_B = B² / (2μ₀)       (magnetic energy density)
```

In an EM wave, u_E = u_B exactly (equipartition).

### Skin Depth

```
δ = √(2 / (ωμσ))
```

where ω = 2πf, μ = μ₀μᵣ, σ = conductivity. Higher frequency → thinner skin depth → more surface confinement.

### Poynting Vector

```
S = (1/μ₀) E × B
```

Direction of energy flow. Time-averaged intensity: `I = E₀² / (2μ₀c)`.

### Boris Algorithm

Given E, B, and time step dt:

```
v⁻ = vⁿ + (q/2m) E dt
t  = (q/2m) B dt
s  = 2t / (1 + |t|²)
v' = v⁻ + v⁻ × t
v⁺ = v⁻ + v' × s
vⁿ⁺¹ = v⁺ + (q/2m) E dt
xⁿ⁺¹ = xⁿ + vⁿ⁺¹ dt
```

The key insight is the rotation step (v⁻ → v⁺): it exactly preserves |v| regardless of dt.

---

## License

MIT
