using LinearAlgebra

# standard triangular lattice vectors
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
            # nearest neighbors (horizontal)
            conditional_add_bond!(bonds, current_idx, to_idx(i + 1, j), standard_spring)
            conditional_add_bond!(bonds, current_idx, to_idx(i - 1, j), standard_spring)
            # nearest neighbors (diagonals)
            if isodd(j) # rows generated with no x-offset
                conditional_add_bond!(bonds, current_idx, to_idx(i - 1, j - 1), standard_spring)
                conditional_add_bond!(bonds, current_idx, to_idx(i,     j - 1), standard_spring)
                conditional_add_bond!(bonds, current_idx, to_idx(i - 1, j + 1), standard_spring)
                conditional_add_bond!(bonds, current_idx, to_idx(i,     j + 1), standard_spring)
            else # iseven(j) -> rows generated with x-offset
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
    WaveguideLattice(nx, ny, kB, Sxx; sliding_dim=0, longitudinal_constraint_dim=0)

Constructs a 2D triangular lattice (lattice constant a = 1, particle mass m = 1)
with harmonic A, B, and C bonds tuned to the constitutively odd state.

# Arguments
- `nx::Int, ny::Int`: Number of unit cells in x and y (`ny` is forced to be even).
- `kB::Real`: Stiffness of the diagonal B bonds.
- `Sxx::Real`: Physical Cauchy prestress (related to dimensionless S by `Sxx = S * kB / √3`).
- `sliding_dim::Int=0`: If 1 or 2, applies sliding waveguide boundaries along that axis.
- `longitudinal_constraint_dim::Int=0`: If 1 or 2, restricts velocity to that axis.
"""
function WaveguideLattice(nx, ny, kB, Sxx; sliding_dim::Int=0, longitudinal_constraint_dim::Int=0)
    if !iseven(ny)
        ny += 1
    end
    kA = -kB / 2
    kC = -kB / 2
    rB = 0.0
    rA = Sxx * sqrt(3.0) / kB
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

            # horizontal (Type A)
            if is_bond_vector(dr, HORIZONTAL_BOND_VEC)
                push!(bonds, Bond(i, j, springA))
            
            # diagonals (Type B)
            elseif is_bond_vector(dr, POSITIVE_DIAGONAL_BOND_VEC) || 
                   is_bond_vector(dr, NEGATIVE_DIAGONAL_BOND_VEC)
                push!(bonds, Bond(i, j, springB))
            
            # vertical NNN (Type C)
            elseif is_bond_vector(dr, VERTICAL_NNN_BOND_VEC)
                push!(bonds, Bond(i, j, springC))
            end
        end
    end

    return system, FixedTopologyCalculator(bonds)
end

"""
    BistableWaveguideLattice(nx, ny, uA, uC, kB, Sxx; sliding_dim=0, longitudinal_constraint_dim=0)

Constructs a triangular waveguide lattice where the negative-stiffness A and C bonds
are replaced by `BistableQuarticSpring` elements with well depths `uA` and `uC`,
matched to the target force and incremental stiffness at the equilibrium bond lengths.
"""
function BistableWaveguideLattice(nx, ny, uA, uC, kB, Sxx; sliding_dim::Int=0, longitudinal_constraint_dim::Int=0)
function BistableWaveguideLattice(nx, ny, uA, uC, kB, Sxx; sliding_dim::Int=0, longitudinal_constraint_dim::Int=0)
    if !iseven(ny)
        ny += 1
    end
    # these are the harmonic properties we're trying to emulate
    kA = -kB / 2
    kC = -kB / 2
    rB = 0.0
    rA_harm = Sxx * sqrt(3.0) / kB
    rC_harm = -rA_harm / sqrt(3.0)

    # relevant lattice distances and target forces
    L_A = 1.0
    L_C = sqrt(3.0)
    
    F_A = kA * (L_A - rA_harm)
    F_C = kC * (L_C - rC_harm)

    (r0A, deltaA) = match_bistable_to_harmonic(L_A, kA, F_A, uA)
    (r0C, deltaC) = match_bistable_to_harmonic(L_C, kC, F_C, uC)

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

            # horizontal (Type A)
            if is_bond_vector(dr, HORIZONTAL_BOND_VEC)
                push!(bonds, Bond(i, j, springA))
            
            # diagonals (Type B)
            elseif is_bond_vector(dr, POSITIVE_DIAGONAL_BOND_VEC) || 
                   is_bond_vector(dr, NEGATIVE_DIAGONAL_BOND_VEC)
                push!(bonds, Bond(i, j, springB))
            
            # vertical NNN (Type C)
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
    imposeSpinWave!(system, modeIndexX, A,  omega)

Modifies positions and velocities to set up an x-propagating circularly polarized wave.
Note this only works for S=0.
"""
function imposeSpinWave!(system::ParticleSystem{D,T,B}, modeIndexX, A, omega) where {D,T,B}
    L = system.boundary_conditions.box_dimensions
    
    qx = 2.0 * pi * modeIndexX / L[1]

    for i in eachindex(system.particles)
        p = system.particles[i]
        phase = qx * p.position[1]
        displacement = SVector{D,T}(A * cos(phase), A * sin(phase))
        new_position = p.position + displacement
        new_velocity = SVector{D,T}(omega * A * sin(phase), -omega * A * cos(phase))
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
