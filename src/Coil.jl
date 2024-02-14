struct Coil
    r::Float64
    α::Float64
    current::Float64
    #buffer_32::Matrix{Float64}
    buffer_32::MMatrix{3, 2}
    buffer_3::MVector{3}
end

function Coil(r, α, current) 
    return Coil(
        r, α, current,
        zeros(MMatrix{3, 2}),
        zeros(MVector{3}),
    )
end

struct CoilConfiguration
    coils::Vector{Coil}
    h::Float64
end

function position_and_derivative!(c::Coil, h, z)
    θ = c.α - h * z
    sin_θ, cos_θ = sincos(θ)
    #c.buffer_32[:, :] = 
    #    [[(c.r * cos_θ) (h * c.r * sin_θ)];
    #    [(c.r * sin_θ) (-h * c.r * cos_θ)];
    #    [z 1]]
    c.buffer_32[1, 1] = c.r * cos_θ
    c.buffer_32[2, 1] = c.r * sin_θ
    c.buffer_32[3, 1] = z
    c.buffer_32[1, 2] = h * c.buffer_32[2, 1]
    c.buffer_32[2, 2] = -h * c.buffer_32[1, 1]
    c.buffer_32[3, 2] = 1.0
    
    #return c.buffer_32
    nothing
end

function compute_B(coil::Coil, h, r_eval, zmax=100.0; rtol=1e-4, atol=1e-10)
    function B_integrand(zp)
        # Simpler version with more allocations:
        #data = γ_and_derivative(coil, h, zp)
        #Δr = r_eval - data[:, 1]
        #temp = normsq(Δr) + (1.0e-30)
        #denominator = temp * sqrt(temp)
        #return cross(data[:, 2], Δr) / denominator

        position_and_derivative!(coil, h, zp)
        Δr = @view coil.buffer_32[:, 1]
        Δr -= r_eval
        temp = normsq(Δr) + (1.0e-30)
        inv_denominator = 1.0 / (temp * sqrt(temp))
        cross!(coil.buffer_3, Δr, coil.buffer_32[:, 2])
        return coil.buffer_3 * inv_denominator

        # Does not work:
        #@. coil.buffer_3 *= inv_denominator
        #return coil.buffer_3
    end
    
    integral1, error_estimate = quadgk(B_integrand, -zmax, 0; rtol=rtol, atol=atol)
    integral2, error_estimate = quadgk(B_integrand, 0, zmax; rtol=rtol, atol=atol)
    return μ0 / (4π) * coil.current * (integral1 + integral2)
end

function compute_B(coil_configuration::CoilConfiguration, r_eval, zmax=5 * 2π / coil_configuration.h)
    return sum(compute_B(coil, coil_configuration.h, r_eval, zmax) for coil in coil_configuration.coils)
end

function compute_ψ(coil::Coil, h, r, α, zmax=100.0; rtol=1e-4, atol=1e-11)
    r0 = coil.r

    function ψ_integrand(zp)
        r_r0_cos = r * r0 * cos(α - coil.α + h * zp)
        numerator = h * h * r_r0_cos + 1
        denominator = sqrt(r * r + r0 * r0 - 2 * r_r0_cos + zp * zp)
        return numerator / denominator
    end

    integral1, error_estimate = quadgk(ψ_integrand, -zmax, 0; rtol=rtol, atol=atol)
    integral2, error_estimate = quadgk(ψ_integrand, 0, zmax; rtol=rtol, atol=atol)
    return -μ0 / (4π) * coil.current * (integral1 + integral2)
end

function compute_ψ(coil_configuration::CoilConfiguration, r, α, zmax=3 * 2π / coil_configuration.h)
    return sum(compute_ψ(coil, coil_configuration.h, r, α, zmax) for coil in coil_configuration.coils)
end

function compute_ψ(coil_configuration::CoilConfiguration, x0::Vector{Float64}, y0::Vector{Float64})
    results = zeros(length(y0), length(x0))
    for jx in eachindex(x0)
        x = x0[jx]
        for jy in eachindex(y0)
            y = y0[jy]
            results[jy, jx] = compute_ψ(coil_configuration, sqrt(x * x + y * y), atan(y, x))
        end
    end
    return results
    #return [compute_ψ(coil_configuration, sqrt(x * x + y * y), atan(y, x)) for x in x0; for y in y0]
end

function ψ_plot(coil_configuration::CoilConfiguration, x0, y0; n_levels=30)
    x0 = collect(x0)
    y0 = collect(y0)
    data = compute_ψ(coil_configuration, x0, y0)
    contour(
        x0,
        y0,
        data,
        levels=n_levels,
        label="", 
        aspect_ratio=:equal,
        c=:rainbow,
    )
    # Plot the coil locations
    for coil in coil_configuration.coils
        x = coil.r * cos(coil.α)
        y = coil.r * sin(coil.α)
        plot!([x], [y],
        color=:black,
        markershape=:xcross,
        label="",
        #markerstrokewidth=5,
    )
    end
    xlabel!("x")
    ylabel!("y")
end

function compute_poincare(coil_configuration::CoilConfiguration, x0, y0, nperiods)
    @assert length(x0) == length(y0)
    results = zeros(length(x0), 2, nperiods + 1)

    function d_position_d_z(position, params, z)
        # Always make z in the first period so we don't need to worry as much
        # about the endpoint of integration for Biot-Savart.
        position_mod = @SVector [position[1], position[2], mod(z, 2π / coil_configuration.h)]
        B = compute_B(coil_configuration, position_mod)
        return B[1:2] / B[3]
    end

    for j in eachindex(x0)
        println("Handling initial condition $j")
        u0 = [x0[j], y0[j]]
        tspan = (0.0, nperiods * 2π / coil_configuration.h)
        prob = ODEProblem(d_position_d_z, u0, tspan)
        #sol = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8)
        sol = solve(prob, saveat=2π / coil_configuration.h, reltol = 1e-4)
        #println("$sol")
        #@show sol
        #@show sol[1]
        #@show sol[:, :]
        results[j, :, 1:nperiods + 1] .= sol[:, 1:nperiods + 1]
    end
    return results
end

function poincare_plot(coil_configuration::CoilConfiguration, x0, y0, nperiods)
    data = compute_poincare(coil_configuration, x0, y0, nperiods)
    scatter!(
        data[:, 1, :]', 
        data[:, 2, :]', 
        label="", 
        markerstrokewidth=0, 
        markersize=2,
        aspect_ratio=:equal,
    )
    # Plot the coil locations
    for coil in coil_configuration.coils
        x = coil.r * cos(coil.α)
        y = coil.r * sin(coil.α)
        plot!([x], [y],
        color=:black,
        markershape=:xcross,
        label="",
        markerstrokewidth=5,
    )
    end
    xlabel!("x")
    ylabel!("y")
end