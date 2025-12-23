@testset "Force Laws" begin

    @testset "HarmonicSpring" begin
        k = 10.0
        r0 = 1.0
        spring = HarmonicSpring(k, r0)

        r_vec = SVector(1.0, 0.0) # Magnitude 1.0
        f = spring_force(spring, r_vec)
        @test norm(f) ≈ 0.0 atol=1e-12

        r_vec_stretch = SVector(1.5, 0.0)
        f_stretch = spring_force(spring, r_vec_stretch)
        @test f_stretch ≈ SVector(5.0, 0.0)

        r_vec_compress = SVector(0.5, 0.0)
        f_compress = spring_force(spring, r_vec_compress)
        @test f_compress ≈ SVector(-5.0, 0.0)
    end

    @testset "BistableQuarticSpring" begin
        r0 = 2.0
        delta = 0.5
        U = 1.0
        spring = BistableQuarticSpring(r0, delta, U)

        f_unstable = spring_force(spring, SVector(2.0, 0.0))
        @test norm(f_unstable) ≈ 0.0 atol=1e-12

        f_stable = spring_force(spring, SVector(2.5, 0.0))
        @test norm(f_stable) ≈ 0.0 atol=1e-12

        r_vec = SVector(2.1, 0.0)
        f_check = spring_force(spring, r_vec)
        
        x = 0.1
        expected_mag = 4.0 * spring.k4 * x^3 - 2.0 * spring.k2 * x
        @test f_check[1] ≈ expected_mag atol=1e-10
    end
end
