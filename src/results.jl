function poincare()
    current = 1.0e6
    coils = [
        Coil(1.0, 0.0, current), 
        Coil(1.0, π, current),
    ]
    h = 1.0
    coil_configuration = CoilConfiguration(coils, h)
    x0 = range(0.025, 0.1, length=4)
    y0 = zeros(length(x0))
    nperiods = 100
    poincare_plot(coil_configuration, x0, y0, nperiods)
end

function ψ_plot_1()
    current = 1.0e6
    coils = [
        Coil(1.0, 0.0, current), 
        Coil(1.0, π, current),
    ]
    h = 1.0
    coil_configuration = CoilConfiguration(coils, h)

    x0 = collect(range(-0.3, 0.3, length=44))
    y0 = collect(range(-0.7, 0.7, length=46))
    ψ_plot(coil_configuration, x0, y0)

    x0 = range(0.025, 0.15, length=4)
    y0 = zeros(length(x0))
    nperiods = 100
    poincare_plot(coil_configuration, x0, y0, nperiods)

end

function h_scan()
    
end