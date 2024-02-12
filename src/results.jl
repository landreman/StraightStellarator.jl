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