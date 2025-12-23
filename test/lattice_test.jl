
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
