using BeatPlanes
using MicroSwimmers
using MicroSwimmersPlots
using GLMakie
using JLD2

@load "results/e2c_sims_final.jld2" e2c_sims
@load "results/e2c_helices_final.jld2" helices

function animate_sim(sim::Excavate2ChlamySimulation, N=100; filename=nothing)
    traj = continue_periodic_trajectory(sim.traj, N)
    animate(traj, generate_e2c(sim), filename=filename, framerate=30)
end

# for convergence checking
function e2c_displacement(N_body=237, Q_body=3996, N_f=94, Q_f=522)
    d = 6*sin(π/7)*sin(π/2)
    p = [π/2, d, -11π/12, -π/2, 0.] 
    e2c = generate_e2c(Excavate2ChlamySimulation("test", N_body=N_body, Q_body=Q_body, N_f=N_f, Q_f=Q_f, β=p[1], d=p[2], Ry_1=p[3], Rx_2=p[4], C=p[5]));
    prob = SwimmingTrajectoryProblem(e2c, t_final=1.0, saveat=1.0, eps=0.15)
    solve_problem!(prob)
    res = check_boundary_conditions(prob.swimming_problem, t=1.0)
    display(viz(prob.swimming_problem))
    @info "" median(res) maximum(res) quantile(res, 0.95) mean(res)
    prob.traj.x[end]
end

function plot_e2c!(parent, e2c::MicroSwimmer; 
        rasterize_body=1,
        bodycolor=Makie.wong_colors()[1], 
        flagellacolor=:forestgreen,
        azimuth=7.1428571796417275, 
        elevation=0.5348572082519532, 
        add_axes=false, 
        show_body=true,
        kwargs...
    )
    ax = Axis3(parent[1,1], 
        aspect=:data, 
        azimuth=azimuth, 
        elevation=elevation,
        protrusions=(0,0,0,0)
        # viewmode=:fit
    )
    if add_axes
        x_ax = [Point3(t, 0., 0.) for t in range(-10,10,10)]
        y_ax = [Point3(0., t, 0.) for t in range(-10,10,10)]
        z_ax = [Point3(0., 0., t) for t in range(-10,10,10)]
        for l in [x_ax, y_ax, z_ax]
            lines!(ax, l, color=:black, linestyle=:dash, linewidth=0.5)
        end
    end
        
    hidedecorations!(ax)
    hidespines!(ax)

    location, orientation = e2c.points.location, e2c.points.orientation
   
    if show_body
        viz!(ax, e2c.body, location, orientation; color=bodycolor, rasterize=rasterize_body)
    end

    f1 = e2c.flagella[1]
    f2 = e2c.flagella[2]
    viz!(ax, f1, location, orientation; color=flagellacolor, kwargs...)
    viz!(ax, f2, location, orientation; color=flagellacolor, kwargs...)
    ax
end

plot_e2c(e2c; kwargs...) = begin
    fig = Figure()
    plot_e2c!(fig[1,1], e2c; kwargs...)
    fig
end

# figures
function e2c_helix_plot!(parent)
    g = parent isa GridLayout ? parent : GridLayout(parent)
    ts = 0.0:0.05:500
   
    ax = Axis3(g[1, 1], 
        aspect=:data, 
        # limits=limits, 
        xlabeloffset = 20,
        ylabeloffset = 20,
        zlabeloffset = 20,
        xlabel=L"x \, (\mu\mathrm{m})",
        ylabel=L"y \, (\mu\mathrm{m})",
        zlabel=L"z \, (\mu\mathrm{m})",
        # left, right, bottom, top
        protrusions=(10, 0, 0, 0),
        # limits=(nothing, 100, nothing, nothing, nothing, nothing)
    )
    hidespines!(ax)
    cmap = :viridis
    cg = cgrad(cmap, length(helices), categorical=true)

    hs = [translate_helix(helices[i], [i*6.0, 0., 0.])(ts) for i in eachindex(helices)]

    proj_x = [[Point3f(300, p[2],p[3]) for p in h] for h in hs]
    proj_y = [[Point3f(p[1], 50., p[3]) for p in h] for h in hs]
    proj_z = [[Point3f(p[1],p[2], -50.) for p in h] for h in hs]

    lines!.([proj_x; proj_y; proj_z], color=(:gray, 0.7))
    # selected = sort(vcat([1, 8, 15, 22, 30], [4, 5, 11, 12, 18, 19, 25, 26, 27]))
    for i in eachindex(helices)
        lines!(ax, hs[i], color=cg[i], linewidth=1.5)
    end

    ax
end

e2c_helix_plot() = begin
    fig = Figure()
    e2c_helix_plot!(fig[1,1])
    fig
end


function plot_helix_measurements()
    fig = Figure()
    titles = [
        "axis velocity",
        "radius", 
        "axis polar angle", 
        "curvature", 
        "axis angular velocity",
        "pitch", 
        "axis azimuthal angle", 
        "torsion"
    ]
    
    axes = [Axis(fig[i,j], xlabel=L"\lambda", title=titles[i*j])  for i in 1:4, j in 1:2]
    # for i in eachindex(axes)
    #     axes[i].ylabel = titles[i]
    # end
    ls = range(0,1,30)
    scatterlines!(axes[1,1], ls[2:end-1], axis_velocity.(helices[2:end-1]))
    scatterlines!(axes[1,2], ls[2:end-1], axis_angular_velocity.(helices[2:end-1]))
    scatterlines!(axes[2,1], ls[2:end-1], radius.(helices[2:end-1]))    
    scatterlines!(axes[2,2], ls[2:end-1], pitch.(helices[2:end-1]))
    scatterlines!(axes[3,1], ls[2:end-1], axis_polar_angle.(helices[2:end-1]))
    scatterlines!(axes[3,2], ls[2:end-1], axis_azimuthal_angle.(helices[2:end-1]))
    scatterlines!(axes[4,1], ls[2:end-1], curvature.(helices[2:end-1]))
    scatterlines!(axes[4,2], ls[2:end-1], torsion.(helices[2:end-1]))

    fig
end

function plot_helix_measurements_minimal!(parent)
    g = parent isa GridLayout ? parent : GridLayout(parent)
    titles = [
        "velocity",
        "curvature", 
        "torsion"
    ]

    ylabels = [
        L"v\;(\mu\mathrm{m}/\text{beat})",
        L"\kappa\;(\mu\mathrm{m}^{-1})",
        L"\tau\;(\mu\mathrm{m}^{-1})",
    ]
    
    axes = [Axis(g[i,1], 
        xlabel= i==3 ? L"\alpha" : "", 
        ylabel = ylabels[i],
        title=titles[i],    
        xticklabelsvisible=i==3 ? true : false
    )  for i in 1:3]

    ls = range(0,1,30)
    ms = 5
    measurements = [
        axis_velocity,
        curvature,
        torsion
    ]
    for i in 1:3
        scatterlines!(axes[i], ls[2:end-1], measurements[i].(helices[2:end-1]),
            linewidth=0.8,
            markersize=ms,
            color=(:navy, 0.7)
        )
    end
    rowgap!(g, 5)

end

plot_helix_measurements_minimal() = begin
    fig = Figure(size=(240,460))
    plot_helix_measurements_minimal!(fig[1,1])
    fig
end

