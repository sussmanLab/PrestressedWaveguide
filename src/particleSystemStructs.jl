""" 
A particle has mass, along with a D-dimensional position and velocity.
"""
struct Particle{D,T}
    position::SVector{D,T}
    velocity::SVector{D,T}
    mass::T
end

function random_particle(box_dimensions::SVector{D,T}, temperature) where {D,T}
    position = rand(SVector{D,T}) .* box_dimensions
    sigma_v = sqrt(T(temperature))
    velocity = randn(SVector{D,T}) .* sigma_v
    mass = one(T)
    return Particle(position, velocity, mass)
end


"""
    Bond(i::Int, j::Int, spring::S)

Represents a single spring connecting particles `i` and `j`.
"""
struct Bond{S}
    i::Int
    j::Int
    spring::S
end

abstract type AbstractBoundaryConditions end

struct OpenBoundaryConditions <: AbstractBoundaryConditions end

struct PeriodicBoundaryConditions{D,T} <:AbstractBoundaryConditions
    box_dimensions::SVector{D,T}
end

struct ParticleSystem{D,T,B}
    particles::Vector{Particle{D,T}}
    boundary_conditions::B
end

"""
    box_displacement(bc, p1_pos, p2_pos)

Calculates the vector pointing from p1 to p2.
For Periodic BCs, this applies the Minimum Image Convention.
"""
function box_displacement(::OpenBoundaryConditions, p1_pos, p2_pos)
    return p2_pos - p1_pos
end

function box_displacement(bc::PeriodicBoundaryConditions{D,T}, p1_pos, p2_pos) where {D,T}
    dr = p2_pos - p1_pos
    
    # We subtract the box dimensions scaled by the nearest integer to dr ./ box_size
    # This effectively wraps the vector components to range [-L/2, L/2]
    return dr .- bc.box_dimensions .* round.(dr ./ bc.box_dimensions)
end

"""
    apply_all_bcs!(system)

Enforces boundary conditions on all particles. 
For periodic systems, this wraps particles back into the box.
"""
function apply_all_bcs!(system::ParticleSystem{D,T,OpenBoundaryConditions}) where {D,T}
    return nothing
end

function apply_all_bcs!(system::ParticleSystem{D,T,PeriodicBoundaryConditions{D,T}}) where {D,T}
    L = system.boundary_conditions.box_dimensions
    for i in eachindex(system.particles)
        p = system.particles[i]
        wrapped_pos = p.position .- L .* floor.(p.position ./ L)
        
        system.particles[i] = Particle(wrapped_pos, p.velocity, p.mass)
    end
end

function kinetic_energy(system::ParticleSystem{D,T,BC}) where {D,T, BC}
    return 0.5 * sum(p.mass * sum(p.velocity.^2) for p in system.particles)
end
