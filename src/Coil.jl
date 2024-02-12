struct Coil
    r::Float64
    α::Float64
    current::Float64
    buffer_32::Matrix{Float64}
end

function Coil(r, α, current) 
    return Coil(
        r, α, current,
        zeros(3, 2),
    )
end

function γ_and_derivative(c::Coil, h, z)
    θ = c.α - h * z
    sin_θ, cos_θ = sincos(θ)
    c.buffer_32[:, :] = 
        [[(c.r * cos_θ) (h * c.r * sin_θ)];
        [(c.r * sin_θ) (-h * c.r * cos_θ)];
        [z 1]]
    
    return c.buffer_32
end

function compute_B(coil::Coil, h, r_eval, zmax=100.0)
    function B_integrand(zp)
        data = γ_and_derivative(coil, h, zp)
        Δr = r_eval - data[:, 1]
        temp = normsq(Δr) + (1.0e-30)
        denominator = temp * sqrt(temp)
        return cross(data[:, 2], Δr) / denominator
    end
    
    integral1, error_estimate = quadgk(B_integrand, -zmax, 0)
    integral2, error_estimate = quadgk(B_integrand, 0, zmax)
    return μ0 / (4π) * coil.current * (integral1 + integral2)
end