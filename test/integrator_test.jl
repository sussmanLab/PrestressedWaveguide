@testset "Integrators" begin

    @testset "Verlet Energy Conservation for oscillator" begin
        p1 = Particle(SVector(0.0, 0.0), SVector(0.0, 0.0), 1.0)
        p2 = Particle(SVector(2.0, 0.0), SVector(0.0, 0.0), 1.0)
        
        box = SVector(10.0, 10.0)
        sys = ParticleSystem([p1, p2], PeriodicBoundaryConditions(box))

        spring = HarmonicSpring(10.0, 1.0)
        
        BondType = Bond{HarmonicSpring{Float64}}
        bonds = [BondType(1, 2, spring)]
        calc = FixedTopologyCalculator(bonds)

        integrator = VerletIntegrator(sys, calc)
        dt = 0.01
        E_initial = total_energy(sys, calc)
        
        for _ in 1:1000
            integrate!(sys, calc, integrator, dt)
        end

        E_final = total_energy(sys, calc)

        @test isapprox(E_initial, E_final, atol=1e-2)
    end
end
