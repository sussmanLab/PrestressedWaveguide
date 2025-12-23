
abstract type AbstractForceCalculator end

"""
    FixedTopologyCalculator(bonds)

A force calculator that iterates over a fixed list of bonds.
Complexity is O(N_bonds).
"""
struct FixedTopologyCalculator{B} <: AbstractForceCalculator
    bonds::Vector{B}
end

"""
    _apply_force_pair!(accelerations, i, j, particles, pairwise_force)

Calculates the force between particles `i` and `j` and updates their
accelerations according to Newton's third law. 
"""
@inline function _apply_force_pair!(accelerations, i, j, system, force)
    p1 = system.particles[i]
    p2 = system.particles[j]

    r_vector = box_displacement(system.boundary_conditions, p1.position, p2.position)

    f_on_p1 = force(r_vector)

    accelerations[i] += f_on_p1 / p1.mass
    accelerations[j] -= f_on_p1 / p2.mass
end

"""
    compute_acceleration!(accelerations, system, calculator::FixedTopologyCalculator)

Iterates only over the defined bonds to compute forces.
"""
function compute_acceleration!(accelerations::Vector{SVector{D,T}}, system::ParticleSystem{D,T,BC}, calculator::FixedTopologyCalculator{B}) where {D,T,BC,B}
    
    fill!(accelerations, zero(SVector{D,T}))

    for bond in calculator.bonds
        _apply_force_pair!(accelerations, bond.i, bond.j, system, r -> spring_force(bond.spring, r))
    end
    
    return nothing
end

function potential_energy(system::ParticleSystem{D,T,BC}, calculator::FixedTopologyCalculator{B}) where {D,T,BC,B}
    pe = 0.0
    for bond in calculator.bonds
        p1 = system.particles[bond.i]
        p2 = system.particles[bond.j]
        dr = box_displacement(system.boundary_conditions, p1.position, p2.position)
        pe += potential_energy(bond.spring, dr)
    end
    return pe
end

function total_energy(system::ParticleSystem{D,T,BC}, calculator::FixedTopologyCalculator{B}) where {D,T,BC,B}
    return kinetic_energy(system) + potential_energy(system, calculator)
end
