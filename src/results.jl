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
    nperiods = 200
    poincare_plot(coil_configuration, x0, y0, nperiods)

end

function ψ_plot_2()
    current = 1.0e6
    coils = [
        Coil(1.0, 0.0, current), 
        Coil(1.0, π, current),
    ]
    h = 1.3
    coil_configuration = CoilConfiguration(coils, h)

    x0 = collect(range(-0.8, 0.8, length=44))
    y0 = collect(range(-1.8, 1.8, length=46))
    fig = ψ_plot(coil_configuration, x0, y0)

    x0 = range(0.1, 0.4, length=4)
    y0 = zeros(length(x0))
    nperiods = 150
    poincare_plot(coil_configuration, x0, y0, nperiods)
    fig
end

function h_scan_2_coils()
    #h_vals = [0.3, 1.0, 3.0]
    h_vals = 10 .^ range(-0.5, 2, length=12)

    current = 1.0e6
    coils = [
        Coil(1.0, 0.0, current), 
        Coil(1.0, π, current),
    ]

    # Set up grid of subplots
    n_plots = length(h_vals)
    n_cols = Int(ceil(0.9 * sqrt(n_plots)))
    n_rows = Int(ceil(n_plots / n_cols))
    @show n_plots, n_rows, n_cols

    layout = (n_rows, n_cols)
    plots = Array{Any, 1}(undef, n_plots)
    scalefontsizes()
    scalefontsizes(0.5)

    Plots.gr_cbar_width[] = 0.005
    
    for jplot in eachindex(h_vals)
        h = h_vals[jplot]
        coil_configuration = CoilConfiguration(coils, h)

        x0 = collect(range(-0.8, 0.8, length=54))
        y0 = collect(range(-0.8, 0.8, length=66))
        plots[jplot] = ψ_plot(coil_configuration, x0, y0; n_levels=50)
        title!("h = $h")
    
    end
    plot(plots..., layout=layout, dpi=100, size=(1100, 850))
    #annotate!((0., 0., directory * filename), subplot=1)
    savefig("20240213-01_h_scan_2_coils.pdf")

end


function h_scan_1_coil()
    #h_vals = [0.3, 1.0, 3.0]
    h_vals = 10 .^ range(-0.5, 1.0, length=12)

    current = 1.0e6
    coils = [
        Coil(1.0, 0.0, current), 
        #Coil(1.0, π, current),
    ]

    # Set up grid of subplots
    n_plots = length(h_vals)
    n_cols = Int(ceil(0.9 * sqrt(n_plots)))
    n_rows = Int(ceil(n_plots / n_cols))
    @show n_plots, n_rows, n_cols

    layout = (n_rows, n_cols)
    plots = Array{Any, 1}(undef, n_plots)
    scalefontsizes()
    scalefontsizes(0.5)

    Plots.gr_cbar_width[] = 0.005
    
    for jplot in eachindex(h_vals)
        h = h_vals[jplot]
        coil_configuration = CoilConfiguration(coils, h)

        xmax = 2.0
        x0 = collect(range(-xmax, xmax, length=54))
        y0 = collect(range(-xmax, xmax, length=66))
        plots[jplot] = ψ_plot(coil_configuration, x0, y0; n_levels=50)
        title!("h = $h")
    
    end
    plot(plots..., layout=layout, dpi=100, size=(1100, 850))
    #annotate!((0., 0., directory * filename), subplot=1)
    savefig("20240213-01_h_scan_1_coil.pdf")

end


function sandbox()
    h = 1.3

    current = 1.0e6
    coils = [
        Coil(1.0, 0.0, current), 
        #Coil(0.5, 0.0, -current), 
        Coil(1.0, π, current),
        #Coil(1.0, π/2, current),
        #Coil(1.0, 3π/2, current),
        #Coil(0.6, π, 0.7*current),
        #Coil(0.5, π, 0.3 * current),
    ]
    coil_configuration = CoilConfiguration(coils, h)

    xmax = 1.0
    ymax = 1.8
    x0 = range(-xmax, xmax, length=154)
    y0 = range(-ymax, ymax, length=166)

    ψ_plot(coil_configuration, x0, y0; n_levels=150)
    plot!(size=(800, 800))
end