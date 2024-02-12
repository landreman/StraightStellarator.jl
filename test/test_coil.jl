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

    @testset "A coil at r=0 should make the expected field for a straight wire" begin
        current = 3.2e6
        for α0 in range(-2.2, 9.6, length=4)
            coil = Coil(0, α0, current)
            for r in range(0.1, 2.6, length=3)
                for z in range(-5, 4, length=5)
                    for θ in range(-2.1, 11.6, length=3)
                        for h in 1:3
                            r_eval = [r * cos(θ), r * sin(θ), z]
                            B = compute_B(coil, h, r_eval)
                            factor = μ0 * current / (2π * r)
                            @test B ≈ [-factor * sin(θ), factor * cos(θ), 0] rtol=1e-3
                        end
                    end
                end
            end
        end
    end

    @testset "Field from a coil list should be independent of the list order" begin
        coil1 = Coil(2.3, 1.3, 1e6)
        coil2 = Coil(1.9, -0.3, 1.2e6)
        h = 1.4
        config1 = CoilConfiguration([coil1, coil2], h)
        config2 = CoilConfiguration([coil2, coil1], h)

        B1 = compute_B(config1, [0.2, 0.1, 0.3])
        B2 = compute_B(config2, [0.2, 0.1, 0.3])
        @test B1 ≈ B2
    end

    @testset "Test computation of Poincare data" begin
        current = 1.0e6
        coils = [
            Coil(1.0, 0.0, current), 
            Coil(1.0, π, current),
        ]
        h = 1.0
        coil_configuration = CoilConfiguration(coils, h)
        x0 = range(0.025, 0.1, length=2)
        y0 = zeros(length(x0))
        nperiods = 20
        compute_poincare(coil_configuration, x0, y0, nperiods)    
    end

    @testset "Test Poincare plot" begin
        current = 1.0e6
        coils = [
            Coil(1.0, 0.0, current), 
            Coil(1.0, π, current),
        ]
        h = 1.0
        coil_configuration = CoilConfiguration(coils, h)
        x0 = range(0.025, 0.1, length=2)
        y0 = zeros(length(x0))
        nperiods = 20
        poincare_plot(coil_configuration, x0, y0, nperiods)    
    end
end
