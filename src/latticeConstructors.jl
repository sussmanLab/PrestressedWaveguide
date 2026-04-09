using LinearAlgebra

# Standard triangular lattice vectors
const A1_vec = SVector{2,Float64}(1.0, 0.0)
const A2_vec = SVector{2,Float64}(0.5, 0.5 * sqrt(3.0))

# checks for constructing the specific waveguide system
const HORIZONTAL_BOND_VEC = A1_vec
const POSITIVE_DIAGONAL_BOND_VEC = A2_vec
const NEGATIVE_DIAGONAL_BOND_VEC = A2_vec - A1_vec
const VERTICAL_NNN_BOND_VEC = 2.0 * A2_vec - A1_vec


"""
    OneDChainLattice(nx)

Creates a 1D chain of `nx+1` particles with Periodic Boundary Conditions.
Returns `(system, calculator)`.
"""
function OneDChainLattice(nx)
    positions = [SVector{2,Float64}(float(i), 0.0) for i in 0:nx]
    N = length(positions)

    box_dims = SVector{2,Float64}(Float64(N), 1.0)
    bc = PeriodicBoundaryConditions(box_dims)

    particles = [Particle(pos, zero(SVector{2,Float64}), 1.0) for pos in positions]
    system = ParticleSystem(particles, bc)

    BondType = Bond{HarmonicSpring{Float64}} 
    bonds = Vector{BondType}()

    k_std = 1.0
    r0_std = 1.0
    spring = HarmonicSpring(k_std, r0_std)

    for i in 1:(N-1)
        push!(bonds, Bond(i, i + 1, spring))
    end
    
    push!(bonds, Bond(N, 1, spring))

    return system, FixedTopologyCalculator(bonds)
end


"""
    SimpleTriangularLattice(nx, ny)

Creates a nearest-neighbor triangular lattice. 
Forces `ny` to be even for consistent staggering.
Returns `(system, calculator)`.
"""
function SimpleTriangularLattice(nx, ny)
    # Ensure even number of rows for periodicity consistency
    if !iseven(ny)
        ny += 1
    end

    y_spacing = 0.5 * sqrt(3.0)
    positions = SVector{2,Float64}[]

    for j in 0:(ny-1) 
        x_offset = isodd(j) ? 0.5 : 0.0
        for i in 0:(nx-1)
            push!(positions, SVector(i + x_offset, j * y_spacing))
        end
    end

    N = length(positions)
    box_dims = SVector{2,Float64}(Float64(nx), ny * y_spacing)
    bc = PeriodicBoundaryConditions(box_dims)

    particles = [Particle(pos, zero(SVector{2,Float64}), 1.0) for pos in positions]
    system = ParticleSystem(particles, bc)

    BondType = Bond{HarmonicSpring{Float64}} 
    bonds = Vector{BondType}()
    standard_spring = HarmonicSpring(1.0, 1.0)

    to_idx(i, j) = (mod1(j, ny) - 1) * nx + mod1(i, nx)

    function conditional_add_bond!(bonds, i_idx, j_idx, spring)
        if i_idx < j_idx
            push!(bonds, Bond(i_idx, j_idx, spring))
        end
    end

    for j in 1:ny
        for i in 1:nx
            current_idx = to_idx(i, j)
            # Nearest Neighbors (Horizontal)
            conditional_add_bond!(bonds, current_idx, to_idx(i + 1, j), standard_spring)
            conditional_add_bond!(bonds, current_idx, to_idx(i - 1, j), standard_spring)
            # Nearest Neighbors (Diagonals)
            if isodd(j) # Rows generated with no x-offset
                conditional_add_bond!(bonds, current_idx, to_idx(i - 1, j - 1), standard_spring)
                conditional_add_bond!(bonds, current_idx, to_idx(i,     j - 1), standard_spring)
                conditional_add_bond!(bonds, current_idx, to_idx(i - 1, j + 1), standard_spring)
                conditional_add_bond!(bonds, current_idx, to_idx(i,     j + 1), standard_spring)
            else # iseven(j) -> Rows generated with x-offset
                conditional_add_bond!(bonds, current_idx, to_idx(i,     j - 1), standard_spring)
                conditional_add_bond!(bonds, current_idx, to_idx(i + 1, j - 1), standard_spring)
                conditional_add_bond!(bonds, current_idx, to_idx(i,     j + 1), standard_spring)
                conditional_add_bond!(bonds, current_idx, to_idx(i + 1, j + 1), standard_spring)
            end
        end
    end

    return system, FixedTopologyCalculator(bonds)
end


"""
    WaveguideLattice(nx, ny, kB, Sxx)

Creates a triangular lattice with specific prestressed bonds (A, B, C) 
(all harmonic springs)
assigned based on spatial orientation.
Returns `(system, calculator)`.
"""
function WaveguideLattice(nx, ny, kB, Sxx; sliding_dim::Int=0, longitudinal_constraint_dim::Int=0)
    if !iseven(ny)
        ny += 1
    end
    kA = -kB / 2
    kC = -kB / 2
    rB = 0.0
    rA = Sxx / kB
    rC = -rA / sqrt(3.0)

    springA = HarmonicSpring(kA, rA)
    springB = HarmonicSpring(kB, rB)
    springC = HarmonicSpring(kC, rC)

    y_spacing = 0.5 * sqrt(3.0)
    positions = SVector{2,Float64}[]
    
    for j in 0:ny-1
        x_offset = isodd(j) ? 0.5 : 0.0
        for i in 0:nx-1
            push!(positions, SVector(i + x_offset, j * y_spacing))
        end
    end
    
    N = length(positions)
    box_dims = SVector{2,Float64}(Float64(nx), ny * y_spacing)
    
    if longitudinal_constraint_dim > 0
        bc = LongitudinalConstraintBC(box_dims, longitudinal_constraint_dim)
    elseif sliding_dim > 0
        p_dim = sliding_dim == 1 ? 2 : 1
        s_dim = sliding_dim
        
        min_coord = minimum(p[s_dim] for p in positions)
        max_coord = maximum(p[s_dim] for p in positions)
        
        side_1_idx = findall(p -> isapprox(p[s_dim], min_coord; atol=1e-5), positions)
        side_2_idx = findall(p -> isapprox(p[s_dim], max_coord; atol=1e-5), positions)
        
        bc = SlidingWaveguideBCs(box_dims, s_dim, p_dim, side_1_idx, side_2_idx, min_coord, max_coord)
    else
        bc = PeriodicBoundaryConditions(box_dims)
    end
    
    particles = [Particle(pos, zero(SVector{2,Float64}), 1.0) for pos in positions]
    system = ParticleSystem(particles, bc)

    BondType = Bond{HarmonicSpring{Float64}} 
    bonds = Vector{BondType}()

    function is_bond_vector(dr, target_vec; atol=0.1)
        return isapprox(dr, target_vec; atol=atol) || isapprox(dr, -target_vec; atol=atol)
    end

    for i in 1:N-1
        for j in i+1:N
            dr = box_displacement(bc, particles[i].position, particles[j].position)

            # Horizontal (Type A)
            if is_bond_vector(dr, HORIZONTAL_BOND_VEC)
                push!(bonds, Bond(i, j, springA))
            
            # Diagonals (Type B)
            elseif is_bond_vector(dr, POSITIVE_DIAGONAL_BOND_VEC) || 
                   is_bond_vector(dr, NEGATIVE_DIAGONAL_BOND_VEC)
                push!(bonds, Bond(i, j, springB))
            
            # Vertical NNN (Type C)
            elseif is_bond_vector(dr, VERTICAL_NNN_BOND_VEC)
                push!(bonds, Bond(i, j, springC))
            end
        end
    end

    return system, FixedTopologyCalculator(bonds)
end

"""
    BistableWaveguideLattice(nx, ny, kB, Sxx)

Creates a triangular lattice with a mixture of harmonic springs and bistable springs, where the bistable springs have been chosen so that their linearization agrees with that of the WaveguideLattice above.
We assume kB = Sxx = 1, and choose the B bond to be a normal harmonic spring. 
We parameterize the A and C bistable springs in terms of the energy well depth of the bistable potential.
Returns `(system, calculator)`.
"""
function BistableWaveguideLattice(nx, ny, uA, uC; sliding_dim::Int=0, longitudinal_constraint_dim::Int=0)
    if !iseven(ny)
        ny += 1
    end
    kB = 1.0
    rB = 0.0

    r0A = 1.0
    deltaA = 2.0 * sqrt(2.0*uA)

    #positive-delta solution for linearization of C bond about sqrt(3)
    function c_bond_helper(u) # bizarre result is an automated common sub-expr reduction
        c225 = u^3
        c226 = 4.0*u
        c227 = -3.0 + c226
        c228 = u*c227
        c229 = 24.0 + c228
        c230 = c225*c229
        c238 = sqrt(c230)
        c222 = u*u
        c223 = 9.0 + u
        c224 = c222*c223
        c239 = 3.0*c238
        c240 = c224 + c239
        c243 = c240^(1/3)
        c249 = c240^(2/3)
        r0=(-14.0*c225 + 2.0*c238*c243 + 2.0*c222*(12.0 + 5.0*c243) - u*(2.0*c238 + 5.0*c249))/(sqrt(3.0)*u*c249)
        delta=  2.0*sqrt(2.0)*sqrt((c222 + c249 - u*(6.0 + c243))/c243)

        return (r0,delta)
    end
    (r0C, deltaC) = c_bond_helper(uC)

    springA = BistableQuarticSpring(r0A, deltaA, uA)
    springB = HarmonicSpring(kB, rB)
    springC = BistableQuarticSpring(r0C, deltaC, uC)

    y_spacing = 0.5 * sqrt(3.0)
    positions = SVector{2,Float64}[]
    
    for j in 0:ny-1
        x_offset = isodd(j) ? 0.5 : 0.0
        for i in 0:nx-1
            push!(positions, SVector(i + x_offset, j * y_spacing))
        end
    end
    
    N = length(positions)
    box_dims = SVector{2,Float64}(Float64(nx), ny * y_spacing)
    
    if longitudinal_constraint_dim > 0
        bc = LongitudinalConstraintBC(box_dims, longitudinal_constraint_dim)
    elseif sliding_dim > 0
        p_dim = sliding_dim == 1 ? 2 : 1
        s_dim = sliding_dim
        
        min_coord = minimum(p[s_dim] for p in positions)
        max_coord = maximum(p[s_dim] for p in positions)
        
        side_1_idx = findall(p -> isapprox(p[s_dim], min_coord; atol=1e-5), positions)
        side_2_idx = findall(p -> isapprox(p[s_dim], max_coord; atol=1e-5), positions)
        
        bc = SlidingWaveguideBCs(box_dims, s_dim, p_dim, side_1_idx, side_2_idx, min_coord, max_coord)
    else
        bc = PeriodicBoundaryConditions(box_dims)
    end
    
    particles = [Particle(pos, zero(SVector{2,Float64}), 1.0) for pos in positions]
    system = ParticleSystem(particles, bc)

    BondType = Union{Bond{HarmonicSpring{Float64}},Bond{BistableQuarticSpring{Float64}}} 
    bonds = Vector{BondType}()

    function is_bond_vector(dr, target_vec; atol=0.1)
        return isapprox(dr, target_vec; atol=atol) || isapprox(dr, -target_vec; atol=atol)
    end

    for i in 1:N-1
        for j in i+1:N
            dr = box_displacement(bc, particles[i].position, particles[j].position)

            # Horizontal (Type A)
            if is_bond_vector(dr, HORIZONTAL_BOND_VEC)
                push!(bonds, Bond(i, j, springA))
            
            # Diagonals (Type B)
            elseif is_bond_vector(dr, POSITIVE_DIAGONAL_BOND_VEC) || 
                   is_bond_vector(dr, NEGATIVE_DIAGONAL_BOND_VEC)
                push!(bonds, Bond(i, j, springB))
            
            # Vertical NNN (Type C)
            elseif is_bond_vector(dr, VERTICAL_NNN_BOND_VEC)
                push!(bonds, Bond(i, j, springC))
            end
        end
    end

    return system, FixedTopologyCalculator(bonds)
end
"""
    imposeWave!(system, modeIndexX, modeIndexY, Ax, Ay)

Modifies the positions of particles to set up a static displacement wave.
Updates `system.particles` in place.
""" 
function imposeWave!(system::ParticleSystem{D,T,B}, modeIndexX, modeIndexY, Ax, Ay) where {D,T,B}
    L = system.boundary_conditions.box_dimensions
    
    qx = 2.0 * pi * modeIndexX / L[1]
    qy = 2.0 * pi * modeIndexY / L[2]

    for i in eachindex(system.particles)
        p = system.particles[i]
        phase = qx * p.position[1] + qy * p.position[2]
        displacement = SVector{D,T}(Ax * cos(phase), Ay * cos(phase))
        new_position = p.position + displacement
        system.particles[i] = Particle(new_position, p.velocity, p.mass)
    end
end

"""
    imposeTravelingWave!(system, modeIndexX, modeIndexY, Ax, Ay, omega)

Modifies positions and velocities to set up a traveling wave.
"""
function imposeTravelingWave!(system::ParticleSystem{D,T,B}, modeIndexX, modeIndexY, Ax, Ay, omega) where {D,T,B}
    L = system.boundary_conditions.box_dimensions
    
    qx = 2.0 * pi * modeIndexX / L[1]
    qy = 2.0 * pi * modeIndexY / L[2]

    for i in eachindex(system.particles)
        p = system.particles[i]
        phase = qx * p.position[1] + qy * p.position[2]
        displacement = SVector{D,T}(Ax * cos(phase), Ay * cos(phase))
        new_position = p.position + displacement
        new_velocity = SVector{D,T}(omega * Ax * sin(phase), omega * Ay * sin(phase))
        system.particles[i] = Particle(new_position, new_velocity, p.mass)
    end
end

"""
    imposeTravelingWavePacket!(system, modeIndexX, modeIndexY, Ax, Ay, omega, sigmaX, sigmaY)

Modifies positions and velocities to set up a Gaussian wave packet centered in the box.
Adds the wave velocity to any existing velocity rather than overwriting it.
"""
function imposeTravelingWavePacket!(system::ParticleSystem{D,T,B}, modeIndexX, modeIndexY, Ax, Ay, omega, sigmaX, sigmaY) where {D,T,B}
    L = system.boundary_conditions.box_dimensions
    
    qx = 2.0 * pi * modeIndexX / L[1]
    qy = 2.0 * pi * modeIndexY / L[2]
    
    centerX = L[1] / 2.0
    centerY = L[2] / 2.0

    for i in eachindex(system.particles)
        p = system.particles[i]
        
        posx = p.position[1]
        posy = p.position[2]
        phase = qx * posx + qy * posy

        gaussianEnvelopeArg = 0.0
        if sigmaX > 0.0
            gaussianEnvelopeArg += ((posx - centerX)^2) / (2.0 * sigmaX^2)
        end
        if sigmaY > 0.0
            gaussianEnvelopeArg += ((posy - centerY)^2) / (2.0 * sigmaY^2)
        end
        gaussianEnvelope = exp(-gaussianEnvelopeArg)

        displacement = gaussianEnvelope * SVector{D,T}(Ax * cos(phase), Ay * cos(phase))
        new_position = p.position + displacement

        wave_velocity = gaussianEnvelope * SVector{D,T}(omega * Ax * sin(phase), omega * Ay * sin(phase))
        new_velocity = p.velocity + wave_velocity

        system.particles[i] = Particle(new_position, new_velocity, p.mass)
    end
end
