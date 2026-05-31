use crate::constants::K_COULOMB;
use crate::electrostatics::{PointCharge, electric_field_point};
use crate::vector::{distance, Vec3};

/// An agent as a charged entity in an EM field analogy.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Agent {
    pub id: String,
    pub position: Vec3,
    pub influence_charge: f64,
    pub sensitivity: f64,
    pub velocity: Vec3,
}

impl Agent {
    pub fn new(id: &str, position: Vec3, influence_charge: f64, sensitivity: f64) -> Self {
        Self { id: id.to_string(), position, influence_charge, sensitivity, velocity: Vec3::zeros() }
    }

    pub fn as_charge(&self) -> PointCharge {
        PointCharge::new(self.position, self.influence_charge)
    }

    pub fn field_at(&self, point: Vec3) -> Vec3 {
        electric_field_point(&self.as_charge(), point)
    }

    pub fn force_from(&self, other: &Agent) -> Vec3 {
        let q1 = PointCharge::new(other.position, other.influence_charge);
        let q2 = PointCharge::new(self.position, self.influence_charge * self.sensitivity);
        crate::electrostatics::coulomb_force(&q1, &q2)
    }

    pub fn total_field(&self, agents: &[Agent]) -> Vec3 {
        let charges: Vec<PointCharge> = agents.iter()
            .filter(|a| a.id != self.id)
            .map(|a| a.as_charge())
            .collect();
        crate::electrostatics::electric_field_superposition(&charges, self.position)
    }

    pub fn apply_field(&mut self, field: Vec3, dt: f64) {
        let force = self.influence_charge * self.sensitivity * field;
        self.velocity += force * dt;
    }

    pub fn step(&mut self, dt: f64) {
        self.position += self.velocity * dt;
    }
}

/// Field map for computing combined influence.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct FieldMap {
    pub agents: Vec<Agent>,
    pub bounds_min: Vec3,
    pub bounds_max: Vec3,
    pub resolution: usize,
}

impl FieldMap {
    pub fn new(agents: Vec<Agent>, bounds_min: Vec3, bounds_max: Vec3, resolution: usize) -> Self {
        Self { agents, bounds_min, bounds_max, resolution }
    }

    pub fn field_at(&self, point: Vec3) -> Vec3 {
        let charges: Vec<PointCharge> = self.agents.iter().map(|a| a.as_charge()).collect();
        crate::electrostatics::electric_field_superposition(&charges, point)
    }

    pub fn potential_at(&self, point: Vec3) -> f64 {
        self.agents.iter().map(|a| {
            let r = distance(a.position, point);
            if r < f64::EPSILON { 0.0 } else { K_COULOMB * a.influence_charge / r }
        }).sum()
    }
}

/// Simulate N-agent interactions.
pub fn simulate_agent_interactions(agents: &mut Vec<Agent>, dt: f64, steps: usize) {
    for _ in 0..steps {
        let forces: Vec<Vec3> = agents.iter().map(|agent| {
            let charges: Vec<PointCharge> = agents.iter()
                .filter(|a| a.id != agent.id)
                .map(|a| a.as_charge())
                .collect();
            let field = crate::electrostatics::electric_field_superposition(&charges, agent.position);
            agent.influence_charge * agent.sensitivity * field
        }).collect();
        for (agent, force) in agents.iter_mut().zip(forces.iter()) {
            agent.velocity += force * dt;
            agent.position += agent.velocity * dt;
        }
    }
}

/// Interaction energy between two agents.
pub fn interaction_energy(a: &Agent, b: &Agent) -> f64 {
    let r = distance(a.position, b.position);
    if r < f64::EPSILON { 0.0 } else { K_COULOMB * a.influence_charge * b.influence_charge / r }
}

/// Total interaction energy.
pub fn total_interaction_energy(agents: &[Agent]) -> f64 {
    let mut energy = 0.0;
    for i in 0..agents.len() {
        for j in (i + 1)..agents.len() {
            energy += interaction_energy(&agents[i], &agents[j]);
        }
    }
    energy
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::vector;

    #[test]
    fn test_agent_field_at() {
        let agent = Agent::new("a1", vector![0.0, 0.0, 0.0], 1e-6, 1.0);
        let field = agent.field_at(vector![1.0, 0.0, 0.0]);
        assert!(field.x > 0.0);
    }

    #[test]
    fn test_agent_force_repulsion() {
        let a1 = Agent::new("a1", vector![0.0, 0.0, 0.0], 1.0, 1.0);
        let a2 = Agent::new("a2", vector![1.0, 0.0, 0.0], 1.0, 1.0);
        assert!(a2.force_from(&a1).x > 0.0);
    }

    #[test]
    fn test_agent_force_attraction() {
        let a1 = Agent::new("a1", vector![0.0, 0.0, 0.0], 1.0, 1.0);
        let a2 = Agent::new("a2", vector![1.0, 0.0, 0.0], -1.0, 1.0);
        assert!(a2.force_from(&a1).x < 0.0);
    }

    #[test]
    fn test_field_map() {
        let agents = vec![Agent::new("a1", vector![0.0, 0.0, 0.0], 1e-6, 1.0)];
        let fm = FieldMap::new(agents, vector![-1.0, -1.0, -1.0], vector![1.0, 1.0, 1.0], 10);
        assert!(fm.field_at(vector![0.5, 0.0, 0.0]).norm() > 0.0);
    }

    #[test]
    fn test_field_map_potential() {
        let agents = vec![Agent::new("a1", vector![0.0, 0.0, 0.0], 1e-6, 1.0)];
        let fm = FieldMap::new(agents, vector![-1.0, -1.0, -1.0], vector![1.0, 1.0, 1.0], 10);
        assert!(fm.potential_at(vector![1.0, 0.0, 0.0]) > 0.0);
    }

    #[test]
    fn test_interaction_energy() {
        let a1 = Agent::new("a1", vector![0.0, 0.0, 0.0], 1.0, 1.0);
        let a2 = Agent::new("a2", vector![1.0, 0.0, 0.0], 1.0, 1.0);
        let expected = K_COULOMB;
        assert!((interaction_energy(&a1, &a2) - expected).abs() / expected < 1e-10);
    }

    #[test]
    fn test_total_interaction_energy() {
        let agents = vec![
            Agent::new("a1", vector![0.0, 0.0, 0.0], 1.0, 1.0),
            Agent::new("a2", vector![1.0, 0.0, 0.0], 1.0, 1.0),
            Agent::new("a3", vector![0.0, 1.0, 0.0], 1.0, 1.0),
        ];
        assert!(total_interaction_energy(&agents) > 0.0);
    }

    #[test]
    fn test_simulate_interactions() {
        let mut agents = vec![
            Agent::new("a1", vector![-1.0, 0.0, 0.0], 1.0, 1.0),
            Agent::new("a2", vector![1.0, 0.0, 0.0], 1.0, 1.0),
        ];
        let d0 = distance(agents[0].position, agents[1].position);
        simulate_agent_interactions(&mut agents, 0.001, 100);
        assert!(distance(agents[0].position, agents[1].position) > d0);
    }

    #[test]
    fn test_agent_superposition() {
        let a1 = Agent::new("a1", vector![1.0, 0.0, 0.0], 1e-6, 1.0);
        let a2 = Agent::new("a2", vector![-1.0, 0.0, 0.0], 1e-6, 1.0);
        let probe = Agent::new("probe", vector![0.0, 0.0, 0.0], 1.0, 1.0);
        let total = probe.total_field(&[a1, a2]);
        assert!(total.x.abs() < 1e-6);
    }
}
