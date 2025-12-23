
@testset "Particle System & Boundaries" begin
    
    @testset "Periodic BC Minimum Image" begin
        box_size = SVector(10.0, 10.0)
        bc = PeriodicBoundaryConditions(box_size)

        p1 = SVector(1.0, 1.0)
        p2 = SVector(9.0, 1.0) 

        dr = box_displacement(bc, p1, p2)
        @test dr ≈ SVector(-2.0, 0.0)

        p3 = SVector(5.0, 1.0)
        p4 = SVector(5.0, 9.5)
        dr_y = box_displacement(bc, p3, p4)
        @test dr_y ≈ SVector(0.0, -1.5)
    end

    @testset "Apply BCs (Wrapping)" begin
        box_size = SVector(10.0, 10.0)
        bc = PeriodicBoundaryConditions(box_size)
        
        p_out = Particle(SVector(11.5, -1.5), SVector(0.0,0.0), 1.0)
        particles = [p_out]
        sys = ParticleSystem(particles, bc)

        apply_all_bcs!(sys)

        @test sys.particles[1].position[1] ≈ 1.5
        @test sys.particles[1].position[2] ≈ 8.5
    end
end
