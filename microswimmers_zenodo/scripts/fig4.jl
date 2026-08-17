include("theme.jl")
include("excavate2chlamy_analysis.jl")
CairoMakie.activate!()


e2c_to_plot = [generate_e2c(e2c_sims[i]) for i in 1:7:28]
push!(e2c_to_plot, generate_e2c(e2c_sims[end]))

function e2c_layout(;draft=false)
    fig = Figure(size = (APS_FULLWIDTH_PT, 350), padding=(0., 0., 0., 0.))

    g = fig[1,1] = GridLayout()

    show_body = !draft

    # first row, morphology
    ga = GridLayout(g[1,1])
    placement = [
        ga[1:2,1:2],
        ga[1,3],
        ga[1,4],
        ga[2,3],
        ga[2,4],
        ga[1:2,5:6]
    ]
    axes = []
    for i in eachindex(e2c_to_plot)
        push!(axes, plot_e2c!(ga[1,i], e2c_to_plot[i]; 
            bodycolor=bc, 
            flagellacolor=fc, 
            show_body=show_body,
            rasterize_body=10
        ))
    end

    rowgap!(ga, 0)
    colgap!(ga, 0)

    Colorbar(g[2,1], 
        limits = (0.0, 1.0), 
        colormap=:viridis, 
        height=5,
        ticks=0:0.25:1,
        label=L"\alpha",
        vertical=false
        # tellwidth=false,
        # tellheight=false
    )

    gb = GridLayout(g[3,1])

    e2c_helix_plot!(gb[1,1])

    plot_helix_measurements_minimal!(gb[1,2])

    rowsize!(g, 1, Relative(0.4))
    colsize!(gb, 2, Relative(0.3))
    rowgap!(g,0)
    colgap!(g,0)
    fig
end
fig = e2c_layout(draft=false)