using StraightStellarator
using Test

@testset "Tests for coil and Biot-Savart functions" begin
    @testset "Biot-Savart should be insensitive to zmax" begin
        coil = Coil(2.3, 1.3, 1e6)

        B1 = compute_B(coil, 1.4, [0.2, 0.1, 0.0])
        B2 = compute_B(coil, 1.4, [0.2, 0.1, 0.0], 1.0e3)
        @show B1
        @show B2
        @test B1 ≈ B2 rtol=1e-3
    end
end
