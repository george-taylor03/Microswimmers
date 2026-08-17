include("theme.jl")
include("chlamy_analysis.jl")
include("filter_feeder_analysis.jl")
include("excavate2chlamy_analysis.jl")
CairoMakie.activate!()


function plot_swimming_prob!(parent, swimmer; kwargs...)
    prob = SwimmingProblem(swimmer, eps=0.15)
    solve_problem!(prob)
    ax = Axis3(parent[1,1], aspect=:data, protrusions=(0,0,0,0))
    hidedecorations!(ax)
    hidespines!(ax)
    viz!(ax, prob; step=4, rasterize=10, kwargs...)
end

@load "results/ff_sim.jld2"
function plot_examples!(parent; kwargs...)
    g = parent isa GridLayout ? parent : GridLayout(parent)

    plot_swimming_prob!(g[1,1], generate_chlamy(chlamy_sims[4,25]); kwargs...)
    Label(g[1,1,Top()], "SwimmingProblem")

    sim = chlamy_sims[8,19]
    chlamy = generate_chlamy(sim)
    traj = continue_periodic_trajectory(sim.traj, 100)
    ax = Axis3(g[1,2], aspect=:data, protrusions=(0,0,0,0))
    Label(g[1,2,Top()], "SwimmingTrajectoryProblem")
    hidedecorations!(ax)
    hidespines!(ax)
    plot_trajectory!(ax, traj, chlamy; kwargs...)

    sim = e2c_sims[5]
    traj = continue_periodic_trajectory(sim.traj, 50)
    ax = Axis3(g[2,2], aspect=:data, protrusions=(0,0,0,0))
    Label(g[2,2,Top()], "Helix")
    hidedecorations!(ax)
    hidespines!(ax)
    lines!(ax, traj.x, color=:forestgreen)
    lines!(ax, helices[5](traj.t), color=:red)

    
    plot_particle_trajectories!(g[2,1], generate_filter_feeder(ff_sim), trajs; kwargs...)
    Label(g[2,1,Top()], "ParticleTrajectoryProblem")
    
    # flow plots
    hm_kwargs=(;colorscale=log10, rasterize=10)
    stream_kwargs=(;linewidth=0.8, arrow_size=4, rasterize=6)
    
    e2c_prob = ResistanceProblem(generate_e2c(e2c_sims[1]), eps=0.15)
    u = FluidVelocity(e2c_prob)
    ax1 = stream!(g[1,3], u, range(-20,20,50), range(-20,20,50); 
        plane=:xz, 
        hm_kwargs=hm_kwargs,
        stream_kwargs=stream_kwargs
    )
    e2c = generate_e2c(e2c_sims[1])
    move_boundary!(e2c, [0. ,0., 0.], rotation_matrix([1., 0., 0.], -π/2), 0.)
    viz!(ax1, e2c; kwargs...)
    Label(g[1,3,Top()], "PlanarVelocityField")
    
    chlamy1 = generate_chlamy(chlamy_sims[1,16])
    chlamy_prob = SwimmingProblem(chlamy1, eps=0.15)
    ave_vf = TimeAveragedPlanarVelocityField(chlamy_prob, range(-30,40,50), range(-30,30,50); num_t=15, period=1.0)
    ax2 = stream!(g[2,3], ave_vf; 
        hm_kwargs=hm_kwargs,
        stream_kwargs=stream_kwargs
    )
    viz!(ax2, chlamy1; kwargs...)
    Label(g[2,3,Top()], "TimeAveragedPlanarVelocityField")
    rowgap!(g, 0)
    colgap!(g, 0)
    colsize!(g, 3, Relative(0.4))
end

function microswimmers_layout(;draft=true)
    fig = Figure(size = (APS_FULLWIDTH_PX, 360), padding=(0., 0., 0., 0.))
    g = fig[1,1] = GridLayout()
    Label(g[1,1,TopLeft()], "A")
    Label(g[1,1,Top()], "microswimmers")
    kwargs = (;rasterize_body=10, bodycolor=bc, color=fc, linewidth=1)
    plot_chlamy!(g[1,2], chlamy_sims[1,16]; kwargs...)
    Label(g[1,2,TopLeft()], "B")
    Label(g[1,2,Top()], "microswimmer models")
    plot_e2c!(g[2,2], generate_e2c(e2c_sims[1]); kwargs...)
    plot_filter_feeder!(g[3,2], ff_sims[8,13]; show_feeding_surface=false, kwargs...)
    plot_examples!(g[1:3,3]; kwargs...)
    Label(g[1:3,3,TopLeft()], "C")
    colsize!(g,3, Relative(0.6))
    colgap!(g, 0)
    fig
end

fig = microswimmers_layout(draft=false)