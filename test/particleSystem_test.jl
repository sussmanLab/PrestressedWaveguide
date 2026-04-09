
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

    @testset "SlidingWaveguideBCs (Y sliding, X periodic)" begin
        box_dims = SVector(10.0, 10.0)
        # 1: X, 2: Y
        bc = SlidingWaveguideBCs(box_dims, 2, 1, [1], [2], 0.0, 10.0)
        
        p1 = SVector(1.0, 1.0)
        p2 = SVector(9.0, 9.0)
        dr = box_displacement(bc, p1, p2)
        @test dr ≈ SVector(-2.0, 8.0)
        
        p_constrained_bottom = Particle(SVector(11.5, -1.0), SVector(1.0, 5.0), 1.0)
        p_constrained_top = Particle(SVector(-1.5, 11.0), SVector(1.0, -5.0), 1.0)
        p_free = Particle(SVector(15.0, 5.0), SVector(2.0, 3.0), 1.0)
        
        sys = ParticleSystem([p_constrained_bottom, p_constrained_top, p_free], bc)
        apply_all_bcs!(sys)
        
        @test sys.particles[1].position[1] ≈ 1.5
        @test sys.particles[2].position[1] ≈ 8.5
        @test sys.particles[3].position[1] ≈ 5.0
        
        @test sys.particles[1].position[2] ≈ 0.0
        @test sys.particles[2].position[2] ≈ 10.0
        @test sys.particles[3].position[2] ≈ 5.0
        
        @test sys.particles[1].velocity[2] ≈ 0.0
        @test sys.particles[2].velocity[2] ≈ 0.0
        @test sys.particles[3].velocity[2] ≈ 3.0
    end

    @testset "SlidingWaveguideBCs (X sliding, Y periodic)" begin
        box_dims = SVector(10.0, 10.0)
        bc = SlidingWaveguideBCs(box_dims, 1, 2, [1], [2], 0.0, 10.0)
        
        p1 = SVector(1.0, 1.0)
        p2 = SVector(9.0, 9.0)
        dr = box_displacement(bc, p1, p2)
        @test dr ≈ SVector(8.0, -2.0)
        
        p_constrained_left = Particle(SVector(-1.0, 11.5), SVector(5.0, 1.0), 1.0)
        p_constrained_right = Particle(SVector(11.0, -1.5), SVector(-5.0, 1.0), 1.0)
        p_free = Particle(SVector(5.0, 15.0), SVector(3.0, 2.0), 1.0)
        
        sys = ParticleSystem([p_constrained_left, p_constrained_right, p_free], bc)
        apply_all_bcs!(sys)
        
        @test sys.particles[1].position[2] ≈ 1.5
        @test sys.particles[2].position[2] ≈ 8.5
        @test sys.particles[3].position[2] ≈ 5.0
        
        @test sys.particles[1].position[1] ≈ 0.0
        @test sys.particles[2].position[1] ≈ 10.0
        @test sys.particles[3].position[1] ≈ 5.0
        
        @test sys.particles[1].velocity[1] ≈ 0.0
        @test sys.particles[2].velocity[1] ≈ 0.0
        @test sys.particles[3].velocity[1] ≈ 3.0
    end

    @testset "LongitudinalConstraintBC" begin
        box_dims = SVector(10.0, 10.0)
        
        bc_x = LongitudinalConstraintBC(box_dims, 1)
        p1 = SVector(1.0, 1.0)
        p2 = SVector(9.0, 9.0)
        dr_x = box_displacement(bc_x, p1, p2)
        @test dr_x ≈ SVector(-2.0, -2.0)
        
        p_test_x = Particle(SVector(11.5, -1.5), SVector(5.0, 3.0), 1.0)
        sys_x = ParticleSystem([p_test_x], bc_x)
        apply_all_bcs!(sys_x)
        @test sys_x.particles[1].position ≈ SVector(1.5, 8.5)
        @test sys_x.particles[1].velocity ≈ SVector(5.0, 0.0)
        
        bc_y = LongitudinalConstraintBC(box_dims, 2)
        dr_y = box_displacement(bc_y, p1, p2)
        @test dr_y ≈ SVector(-2.0, -2.0)
        
        p_test_y = Particle(SVector(11.5, -1.5), SVector(5.0, 3.0), 1.0)
        sys_y = ParticleSystem([p_test_y], bc_y)
        apply_all_bcs!(sys_y)
        @test sys_y.particles[1].position ≈ SVector(1.5, 8.5)
        @test sys_y.particles[1].velocity ≈ SVector(0.0, 3.0)
    end
end
