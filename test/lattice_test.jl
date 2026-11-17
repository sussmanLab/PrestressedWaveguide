
@testset "Lattice Constructors" begin
    
    @testset "OneD Chain" begin
        nx = 5
        sys, calc = OneDChainLattice(nx)
        
        @test length(sys.particles) == nx + 1
        
        @test length(calc.bonds) == nx + 1
    end

    @testset "Triangular Lattice Dimensions" begin
        nx = 4
        ny = 4
        sys, calc = SimpleTriangularLattice(nx, ny)

        @test length(sys.particles) == nx * ny
        
        @test length(calc.bonds) == 3 * length(sys.particles)
    end

    @testset "Waveguide Lattice Parameterized Types" begin
        nx, ny = 4, 4
        kB, Sxx = 10.0, 0.1
        sys, calc = WaveguideLattice(nx, ny, kB, Sxx)

        @test eltype(calc.bonds) <: Bond{<:HarmonicSpring}
        @test isconcretetype(eltype(calc.bonds))
    end
end

@testset "Wave Initialization" begin

    L = 10.0
    box_dims = SVector(L, L)
    bc = PeriodicBoundaryConditions(box_dims)
    
    particles = Vector{Particle{2,Float64}}()
    for y in 0:9
        for x in 0:9
            pos = SVector(float(x), float(y))
            vel = SVector(0.0, 0.0) # Start at rest
            push!(particles, Particle(pos, vel, 1.0))
        end
    end
    
    system = ParticleSystem(particles, bc)

    modeX, modeY = 1.0, 0.0
    Ax, Ay = 0.1, 0.0     
    omega = 1.0
    sigma = 1.0          
    
    imposeTravelingWavePacket!(system, modeX, modeY, Ax, Ay, omega, sigma, sigma)

    p_center = system.particles[56] 
    
    @test p_center.position[1] ≈ 4.9 atol=1e-5
    @test p_center.position[2] ≈ 5.0 atol=1e-5
    
    @test norm(p_center.velocity) < 1e-10
    
    p_edge = system.particles[1]
    
    @test p_edge.position[1] ≈ 0.0 atol=1e-10
    @test p_edge.position[2] ≈ 0.0 atol=1e-10
    @test norm(p_edge.velocity) < 1e-10

    p_q = system.particles[53]
    expected_disp_x = exp(-4.5) * Ax * cos(0.4*pi)
    actual_disp_x = p_q.position[1] - 2.0
    @test actual_disp_x ≈ expected_disp_x atol=1e-5
end
