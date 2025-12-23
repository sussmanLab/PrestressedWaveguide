abstract type AbstractIntegrator end

struct ForwardEuler{D,T} <: AbstractIntegrator 
    accelerations::Vector{SVector{D,T}}
end
function ForwardEuler(system::ParticleSystem{D,T,B}) where {D,T,B}
    return ForwardEuler(zeros(SVector{D,T},length(system.particles)))
end

struct SymplecticEuler{D,T}<: AbstractIntegrator
    accelerations::Vector{SVector{D,T}}
end
function SymplecticEuler(system::ParticleSystem{D,T,B}) where {D,T,B}
    return SymplecticEuler(zeros(SVector{D,T},length(system.particles)))
end

struct VerletIntegrator{D,T} <: AbstractIntegrator
    accelerations::Vector{SVector{D,T}}
end
function VerletIntegrator(system::ParticleSystem{D,T,B},force_calculator) where {D,T,B}
    initial_accelerations = zeros(SVector{D,T}, length(system.particles))
    compute_acceleration!(initial_accelerations, system, force_calculator)
    return VerletIntegrator(initial_accelerations)
end


function euler_step(p::Particle, acceleration, dt)
    new_position = p.position + p.velocity * dt
    new_velocity = p.velocity + acceleration * dt
    return Particle(new_position, new_velocity, p.mass)
end

"""
    integrate!(system, force_calculator, integrator::ForwardEuler, dt)

Advances a system of particles by one `dt` using the Forward Euler integration method.
"""
function integrate!(system::ParticleSystem{D,T,B}, force_calc,
                integrator::ForwardEuler{D,T}, dt::Float64) where {D,T,B}

    compute_acceleration!(integrator.accelerations,system,force_calc)
    system.particles .= euler_step.(system.particles,
                                    integrator.accelerations, dt)
    apply_all_bcs!(system)
end


function symplectic_euler_step(p::Particle, acceleration, dt)
    new_velocity = p.velocity + acceleration * dt
    new_position = p.position + new_velocity * dt
    return Particle(new_position, new_velocity, p.mass)
end
"""
    integrate!(system, force_calculator, integrator::SymplecticEuler, dt)

Advances a system by one `dt` using the Symplectic Euler integration method. The difference here is that we *first* update the velocity, and then update the position according to v(t+Δt), rather than updating the position according to v(t).

("Symplectic Euler" is also sometimes known as "semi-implicit Euler")
"""
function integrate!(system::ParticleSystem{D,T,B}, force_calc,
              integrator::SymplecticEuler{D,T}, dt::Float64) where {D,T,B}

    compute_acceleration!(integrator.accelerations,system,force_calc)
    system.particles .= symplectic_euler_step.(system.particles,
                                    integrator.accelerations, dt)
    apply_all_bcs!(system)
end


"""
    integrate!(system, force_calculator, integrator::VerletIntegrator, dt)

Advances the system by one `dt` using the Velocity Verlet integration method. This is a second-order symplectic integrator known for its high accuracy and excellent long-term energy conservation.

The algorithmic pattern is:
    1. r(t+Δt) = r(t) + Δt*v(t) + 0.5*Δt^2*a(t)
    2. calculate a(t+Δt) using the new positions
    3. v(t+Δt) = v(t) + 0.5*Δt*(a(t) + a(t+Δt))
"""
function integrate!(system::ParticleSystem{D,T,B}, force_calculator,
            integrator::VerletIntegrator{D,T}, dt::Float64) where {D,T,B}

    n = length(system.particles)
    old_accelerations = integrator.accelerations
    # 1. r(t+Δt) = r(t) + Δt*v(t) + 0.5*Δt^2*a(t)
    for i in eachindex(system.particles)
        p = system.particles[i]
        new_position = p.position + p.velocity * dt + 0.5 * dt^2 * old_accelerations[i]
        v_half_step = p.velocity + 0.5 * dt * old_accelerations[i]
        system.particles[i] = Particle(new_position, v_half_step, p.mass)
    end

    apply_all_bcs!(system)

    # 2. calculate a(t+Δt) using the new positions
    compute_acceleration!(integrator.accelerations,system,force_calculator)

    # 3. v(t+Δt) = v(t) + 0.5*Δt*(a(t) + a(t+Δt))
    for i in eachindex(system.particles)
        p = system.particles[i]

        # p.velocity is currently the half-step velocity
        final_velocity = p.velocity + 0.5 * dt * integrator.accelerations[i]

        # Update the particle with its final velocity for this step
        system.particles[i] = Particle(p.position, final_velocity, p.mass)
    end
end

