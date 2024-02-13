using StraightStellarator
using Test

@testset "Tests for utility functions" begin
    @testset "Test cross 2 ways" begin 
        v1 = [0.2, -0.5, -0.3]
        v2 = [-0.1, 1.1, -0.9]
        v3 = cross(v1, v2)
        v4 = zeros(3)
        cross!(v4, v1, v2)
        @test v3 ≈ v4
    end
end
