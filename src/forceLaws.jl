function norm_squared(r::SVector{D,T}) where {D,T}
    return sum(r .* r)
end

"""
    AbstractSpring

An abstract supertype for all spring-like interactions.
"""
abstract type AbstractSpring end

"""
    HarmonicSpring(k::Float64, r0::Float64)

A standard harmonic spring with potential energy `V(r) = 1/2 * k * (r - r0)^2`.
"""
struct HarmonicSpring{T} <: AbstractSpring
    k::T    # Spring constant
    r0::T   # Rest length
end

"""
    BistableQuarticSpring(r0, delta, U)

A double-well spring -- has two minima of depth `U` located at `r0 ± delta`

# Arguments
- `r0::Float64`: The distance corresponding to the local maximum of the potential
- `delta::Float64`: The separation of each minimum from `r0`.
- `U::Float64`: The depth of the potential wells relative to the peak at `r0`.
"""
struct BistableQuarticSpring{T} <: AbstractSpring
    r0::T
    # delta::Float64
    # U::Float64
    
    # Pre-calculated coefficients V(x) = k4*x⁴ - k2*x² where x = r - r0
    k2::T
    k4::T

    function BistableQuarticSpring(r0::T, delta::T, U::T) where {T}
        @assert r0 > 0 && delta > 0 && U > 0 "All parameters must be positive."
        @assert delta < r0 "Separation delta must be smaller than the central distance r0 to ensure the first minimum is at a positive distance."
        
        # Derive internal coefficients from the intuitive parameters
        k2 = T(2.0) * U / (delta^2)
        k4 = U / (delta^4)
        
        new{T}(r0, k2, k4)
    end
end

"""
    spring_force(spring, r_vector)

Calculates the force vector exerted on particle 1 by particle 2,
given the displacement vector `r_vector` pointing from p1 to p2.
"""
function spring_force(s::HarmonicSpring, r_vector)
    r_dist = sqrt(norm_squared(r_vector))
    if r_dist < 1e-12
        return zero(r_vector)
    end
    magnitude = s.k * (r_dist - s.r0)
    return (magnitude / r_dist) * r_vector
end
function spring_force(s::BistableQuarticSpring, r_vector)
    r_dist = sqrt(norm_squared(r_vector))
    if r_dist < 1e-12
        return zero(r_vector)
    end
    x = (r_dist - s.r0)
    force_magnitude = 4.0 * s.k4 * x^3 - 2.0 * s.k2 * x
    return (force_magnitude / r_dist) * r_vector
end

function potential_energy(s::HarmonicSpring, r_vector)
    r_dist = sqrt(norm_squared(r_vector))
    dr = (r_dist - s.r0)
    return 0.5*s.k*dr*dr
end

function potential_energy(s::BistableQuarticSpring, r_vector)
    r_dist = sqrt(norm_squared(r_vector))
    dr = (r_dist - s.r0)
    return s.k4*dr*dr*dr*dr - s.k2*dr*dr
end
