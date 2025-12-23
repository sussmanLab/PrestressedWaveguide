"""
    average_energy_density(system, calculator)

Calculates the total energy  per unit area.
"""
function average_energy_density(system::ParticleSystem{D,T,BC}, calculator) where {D, T,BC}
    volume = prod(system.boundary_conditions.box_dimensions)
    
    return total_energy(system,calculator) / volume
end

"""
    calculate_energy_flux(system, calculator, axis::Int, line_position::Float64)

Calculates the instantaneous energy flux across a line perpendicular to to the `axis` (1 for x-flux, 2 for y, etc), at `line_position` 
"""
function calculate_energy_flux(system::ParticleSystem{2,T,BC}, calculator, axis::Int, line_position::Float64) where {T,BC}
    total_flux = 0.0
    bc = system.boundary_conditions

    for bond in calculator.bonds
        p1 = system.particles[bond.i]
        p2 = system.particles[bond.j]

        dr = box_displacement(bc, p1.position, p2.position)
        
        p1_coord = p1.position[axis]
        p2_virtual_coord = p1_coord + dr[axis]

        if (p1_coord - line_position) * (p2_virtual_coord - line_position) < 0
            f_on_p1 = spring_force(bond.spring, dr)
            contribution = 0.5 * dot(f_on_p1, (p1.velocity + p2.velocity))
            if dr[axis] < 0
                total_flux += contribution
            else
                total_flux -= contribution
            end
        end
    end
    return total_flux
end

"""
    PseudomomentumCalculator

A helper struct to cache the inverse reference matrices for deformation gradient calculations.
This avoids re-computing the reference configuration properties at every time step.
"""
struct PseudomomentumCalculator{D,T}
    # B_inv[i] stores (Σ ΔX_ij ⊗ ΔX_ij)⁻¹ for particle i
    inverse_reference_matrices::Vector{SMatrix{D,D,T,4}} 
    adjacency_list::Vector{Vector{Int}}
    reference_positions::Vector{SVector{D,T}}
    box_size::SVector{D,T}
end

"""
    PseudomomentumCalculator(system, calculator)

Constructs the calculator by analyzing the initial (reference) topology.
"""
function PseudomomentumCalculator(system::ParticleSystem{D,T,BC}, calculator::FixedTopologyCalculator) where {D,T,BC}
    N = length(system.particles)
    
    adj = [Int[] for _ in 1:N]
    for bond in calculator.bonds
        push!(adj[bond.i], bond.j)
        push!(adj[bond.j], bond.i)
    end

    ref_pos = [p.position for p in system.particles]
    inv_mats = Vector{SMatrix{D,D,T,4}}(undef, N)
    
    bc = system.boundary_conditions

    for i in 1:N
        B = zeros(MMatrix{D,D,T})
        for j in adj[i]
            dX = box_displacement(bc, ref_pos[i], ref_pos[j])
            B += dX * dX'
        end
        
        detB = det(B)
        if abs(detB) < 1e-12
            inv_mats[i] = SMatrix{D,D,T}(I)
        else
            inv_mats[i] = SMatrix{D,D,T}(inv(B))
        end
    end

    return PseudomomentumCalculator(inv_mats, adj, ref_pos, bc.box_dimensions)
end

"""
    calculate_pseudomomentum_density(system, pm_calc)

Computes the pseudomomentum density using the precomputed reference matrices.
"""
function calculate_pseudomomentum_density(system::ParticleSystem{D,T,BC}, calc::PseudomomentumCalculator{D,T}) where {D,T,BC}
    total_pseudo = zeros(SVector{D,T})
    bc = system.boundary_conditions

    zero_mat = zero(SMatrix{D,D,T,4})

    for i in eachindex(system.particles)
        p = system.particles[i]
        
        # Calculate A = Σ Δr_ij ⊗ ΔX_ij
        A = zeros(MMatrix{D,D,T})
        
        for j in calc.adjacency_list[i]
            dX = box_displacement(bc, calc.reference_positions[i], calc.reference_positions[j])
            dr = box_displacement(bc, system.particles[i].position, system.particles[j].position)
            A += dr * dX'
        end

        # Deformation Gradient F_i = A * B_inv
        Fi = SMatrix(A) * calc.inverse_reference_matrices[i]

        pi = p.mass * p.velocity

        # Add F_i^T * p_i
        total_pseudo += Fi' * pi
    end
    
    area = calc.box_size[1] * calc.box_size[2]
    return total_pseudo / area
end


"""
    measure_wave_properties(system, calculator, integrator, params...)

Runs the measurement simulation.
"""
function measure_wave_properties(
    system::ParticleSystem{D,T,BC}, 
    force_calculator, 
    integrator,
    modeIndexX, modeIndexY, 
    Ax, Ay, omega,
    tMeasurement, dt; 
    tSettle=20.0
) where {D,T,BC}

    pm_calc = PseudomomentumCalculator(system, force_calculator)
    ground_state_energy_density = average_energy_density(system, force_calculator)
    
    imposeTravelingWave!(system, modeIndexX, modeIndexY, Ax, Ay, omega)
    
    compute_acceleration!(integrator.accelerations, system, force_calculator)

    if tSettle > 0
        run_simulation!(system, force_calculator, integrator, dt, round(Int, tSettle/dt))
    end

    steps = round(Int, tMeasurement / dt)
    if steps == 0 error("tMeasurement too short") end

    flux_sum_x = 0.0
    flux_sum_y = 0.0
    energy_dens_sum = 0.0
    mom_dens_sum = zero(SVector{D,T})

    L = system.boundary_conditions.box_dimensions
    x_line = L[1] / 2.0 + 0.25
    y_line = L[2] / 2.0 + 0.25

    for _ in 1:steps
        integrate!(system, force_calculator, integrator, dt)
        
        flux_sum_x += calculate_energy_flux(system, force_calculator, 1, x_line)
        flux_sum_y += calculate_energy_flux(system, force_calculator, 2, y_line)
        
        curr_e_dens = average_energy_density(system, force_calculator)
        energy_dens_sum += (curr_e_dens - ground_state_energy_density)
        
        mom_dens_sum += calculate_pseudomomentum_density(system, pm_calc)
    end

    avg_flux_x = flux_sum_x / steps
    avg_flux_y = flux_sum_y / steps
    
    jEx = avg_flux_x / L[2] 
    jEy = avg_flux_y / L[1] 
    energy_current_density = SVector(jEx, jEy)
    
    energy_density = energy_dens_sum / steps
    momentum_density = mom_dens_sum / steps
    
    group_velocity = abs(energy_density) > 1e-9 ? energy_current_density / energy_density : zero(SVector{D,T})

    return (
        energy_current_density = energy_current_density,
        energy_density = energy_density,
        momentum_density = momentum_density,
        group_velocity = group_velocity
    )
end
