include("theme.jl")
include("filter_feeder_analysis.jl")
CairoMakie.activate!()


function filter_feeder_layout(; draft=false)
    fig = Figure(size = (APS_FULLWIDTH_PT, 350), padding=(0., 0., 0., 0.))

    g = fig[1,1] = GridLayout()
    g1 = g[1,1] = GridLayout()

    kwargs=(;bodycolor=bc, color=fc, linewidth=1, rasterize_body=10)
    plot_filter_feeder!(g1[1,1], ff_sims[8,16]; kwargs...)
    plot_filter_feeder!(g1[2,1], inverted_sims[7,2]; kwargs...)
    
    g2 = g[1,2] = GridLayout()
    flow_plot_comparison!(g2)
    
    g3 = g[1,3] = GridLayout()
    efficiency_plot!(g3)

    colgap!(g, 0)
    fig
end
fig = filter_feeder_layout(draft=false)