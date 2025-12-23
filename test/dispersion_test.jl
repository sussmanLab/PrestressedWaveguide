function analytical_group_velocity_L(k)
    omega = sqrt(3.0 - 2.0*cos(k) - cos(k/2))
    vg = (2.0*sin(k) + 0.5*sin(k/2)) / (2 * omega)
    return vg
end

function analytical_group_velocity_T(k)
    omega = sqrt(3.0 - 3.0*cos(k/2))
    vg = (1.5 * sin(k/2)) / (2 * omega)
    return vg
end

@testset "Triangular lattice dispersion" begin
    nx = 50
    ny = 8
    
    dt = 0.01
    amplitude = 0.01
    tMeasMultiples = 4.0
    
    modeIdx = 4
    
    sys_ref, _ = SimpleTriangularLattice(nx, ny)
    L_box = sys_ref.boundary_conditions.box_dimensions[1]
    kx = 2 * pi * modeIdx / L_box

    @testset "Longitudinal Mode (L)" begin
        sys, calc = SimpleTriangularLattice(nx, ny)
        integrator = VerletIntegrator(sys, calc)
        omega_L = sqrt(3.0 - cos(kx) - 2.0*cos(kx/2))
        period_L = 2 * pi / omega_L
        tMeas_L = tMeasMultiples * period_L
        
        expected_vg = analytical_group_velocity_L(kx)

        results = measure_wave_properties(
            sys, calc, integrator,
            float(modeIdx), 0.0,
            amplitude, 0.0,    
            omega_L,
            tMeas_L, dt;
            tSettle = 0.0 * dt
        )

        measured_vg_x = results.group_velocity[1]

        println("  L-Mode k=$(round(kx,digits=2)): Expected vg ≈ $(round(expected_vg,digits=4)), Got $(round(measured_vg_x,digits=4))")

        ratio_vg = results.energy_current_density[1] / results.energy_density
        @test isapprox(ratio_vg, measured_vg_x, atol=1e-3)

        @test isapprox(measured_vg_x, expected_vg, rtol=0.05)
    end

    @testset "Transverse Mode (T)" begin
        sys, calc = SimpleTriangularLattice(nx, ny)
        integrator = VerletIntegrator(sys, calc)
        omega_T = sqrt(3.0 - 3.0*cos(kx/2))
        period_T = 2 * pi / omega_T
        tMeas_T = tMeasMultiples * period_T
        
        expected_vg = analytical_group_velocity_T(kx)

        results = measure_wave_properties(
            sys, calc, integrator,
            float(modeIdx), 0.0,
            0.0, amplitude,
            omega_T,
            tMeas_T, dt;
            tSettle = 0.0 * dt
        )

        measured_vg_x = results.group_velocity[1]

        println("  T-Mode k=$(round(kx,digits=2)): Expected vg ≈ $(round(expected_vg,digits=4)), Got $(round(measured_vg_x,digits=4))")

        ratio_vg = results.energy_current_density[1] / results.energy_density
        @test isapprox(ratio_vg, measured_vg_x, atol=1e-3)

        @test isapprox(measured_vg_x, expected_vg, rtol=0.05)
    end
end
