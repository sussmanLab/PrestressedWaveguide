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
function WaveguideLattice(nx, ny, kB, Sxx)
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
    bc = PeriodicBoundaryConditions(box_dims)
    
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

            # Check Horizontal (Type A)
            if is_bond_vector(dr, HORIZONTAL_BOND_VEC)
                push!(bonds, Bond(i, j, springA))
            
            # Check Diagonals (Type B)
            elseif is_bond_vector(dr, POSITIVE_DIAGONAL_BOND_VEC) || 
                   is_bond_vector(dr, NEGATIVE_DIAGONAL_BOND_VEC)
                push!(bonds, Bond(i, j, springB))
            
            # Check Vertical NNN (Type C)
            elseif is_bond_vector(dr, VERTICAL_NNN_BOND_VEC)
                push!(bonds, Bond(i, j, springC))
            end
        end
    end

    return system, FixedTopologyCalculator(bonds)
end
