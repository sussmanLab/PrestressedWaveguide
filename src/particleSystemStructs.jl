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

struct LongitudinalConstraintBC{D,T} <: AbstractBoundaryConditions
    box_dimensions::SVector{D,T}
    longitudinal_dim::Int
end

struct SlidingWaveguideBCs{T} <: AbstractBoundaryConditions
    box_dimensions::SVector{2,T}
    sliding_dim::Int # 1 for x, 2 for y
    periodic_dim::Int # 2 for y, 1 for x
    side_1_indices::Vector{Int}
    side_2_indices::Vector{Int}
    side_1_coord::T
    side_2_coord::T
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
    return dr .- bc.box_dimensions .* round.(dr ./ bc.box_dimensions)
end

function box_displacement(bc::LongitudinalConstraintBC{D,T}, p1_pos, p2_pos) where {D,T}
    dr = p2_pos - p1_pos
    return dr .- bc.box_dimensions .* round.(dr ./ bc.box_dimensions)
end

function box_displacement(bc::SlidingWaveguideBCs{T}, p1_pos, p2_pos) where {T}
    dr = p2_pos - p1_pos
    p_dim = bc.periodic_dim
    dp = dr[p_dim] - bc.box_dimensions[p_dim] * round(dr[p_dim] / bc.box_dimensions[p_dim])
    
    if p_dim == 1
        return SVector{2,T}(dp, dr[2])
    else
        return SVector{2,T}(dr[1], dp)
    end
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

function apply_all_bcs!(system::ParticleSystem{2,T,LongitudinalConstraintBC{2,T}}) where {T}
    L = system.boundary_conditions.box_dimensions
    l_dim = system.boundary_conditions.longitudinal_dim
    for i in eachindex(system.particles)
        p = system.particles[i]
        wrapped_pos = p.position .- L .* floor.(p.position ./ L)
        
        if l_dim == 1
            new_vel = SVector{2,T}(p.velocity[1], zero(T))
        else
            new_vel = SVector{2,T}(zero(T), p.velocity[2])
        end
        
        system.particles[i] = Particle(wrapped_pos, new_vel, p.mass)
    end
end

function apply_all_bcs!(system::ParticleSystem{2,T,SlidingWaveguideBCs{T}}) where {T}
    bc = system.boundary_conditions
    L = bc.box_dimensions
    p_dim = bc.periodic_dim
    s_dim = bc.sliding_dim
    
    for i in eachindex(system.particles)
        p = system.particles[i]
        
        wrapped_p_pos = p.position[p_dim] - L[p_dim] * floor(p.position[p_dim] / L[p_dim])
        
        new_pos_p = wrapped_p_pos
        new_pos_s = p.position[s_dim]
        new_vel_s = p.velocity[s_dim]

        if i in bc.side_1_indices
            new_pos_s = bc.side_1_coord
            new_vel_s = zero(T)
        elseif i in bc.side_2_indices
            new_pos_s = bc.side_2_coord
            new_vel_s = zero(T)
        end
        
        if p_dim == 1
            system.particles[i] = Particle(
                SVector{2,T}(new_pos_p, new_pos_s), 
                SVector{2,T}(p.velocity[1], new_vel_s), 
                p.mass
            )
        else
            system.particles[i] = Particle(
                SVector{2,T}(new_pos_s, new_pos_p), 
                SVector{2,T}(new_vel_s, p.velocity[2]), 
                p.mass
            )
        end
    end
end

function kinetic_energy(system::ParticleSystem{D,T,BC}) where {D,T, BC}
    return 0.5 * sum(p.mass * sum(p.velocity.^2) for p in system.particles)
end
