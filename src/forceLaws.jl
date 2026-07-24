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

"""
    match_bistable_to_harmonic(L, K, F, u)

Solve the cubic equation to find (r0, delta) that match a harmonic spring's
stiffness (K) and force (F) at a specific lattice distance (L). The well depth,
u, is a free parameter.
"""
function match_bistable_to_harmonic(L, K, F, u)
    if abs(F) < 1e-12
        return L, sqrt(-4.0 * u / K)
    end

    a = K^2
    b = -6.0 * K * F
    c = 9.0 * F^2 - 8.0 * u * K
    d_coef = 8.0 * u * F

    p = (3.0 * a * c - b^2) / (3.0 * a^2)
    q = (2.0 * b^3 - 9.0 * a * b * c + 27.0 * a^2 * d_coef) / (27.0 * a^3)

    Delta = (q / 2.0)^2 + (p / 3.0)^3

    xs = Float64[]
    if Delta <= 0
        # use trig if there are three real roots 
        r = sqrt(-(p / 3.0)^3)
        theta = acos(-q / (2.0 * r))
        push!(xs, 2.0 * cbrt(r) * cos(theta / 3.0) - b / (3.0 * a))
        push!(xs, 2.0 * cbrt(r) * cos((theta + 2.0 * pi) / 3.0) - b / (3.0 * a))
        push!(xs, 2.0 * cbrt(r) * cos((theta + 4.0 * pi) / 3.0) - b / (3.0 * a))
    else
        # otherwise Cardano gives us ne real root
        u_p = cbrt(-q / 2.0 + sqrt(Delta))
        v_p = cbrt(-q / 2.0 - sqrt(Delta))
        push!(xs, u_p + v_p - b / (3.0 * a))
    end

    valid_roots = Float64[]
    for x in xs
        denom = K * x - 3.0 * F
        if abs(denom) > 1e-12
            d2 = 8.0 * u * x / denom
            if d2 > 0
                push!(valid_roots, x)
            end
        end
    end

    if isempty(valid_roots)
        error("No valid matching bistable spring parameters found for K=$K, F=$F, u=$u")
    end

    positive_roots = filter(x -> x > 0, valid_roots)
    x_sol = isempty(positive_roots) ? minimum(valid_roots) : minimum(positive_roots)

    r0 = L - x_sol
    delta = sqrt(8.0 * u * x_sol / (K * x_sol - 3.0 * F))

    return r0, delta
end
