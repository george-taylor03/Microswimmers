using BeatPlanes
using MicroSwimmers
using MicroSwimmersPlots
using LaTeXStrings
using JLD2
using GLMakie
using LinearAlgebra

@load "results/chlamy_sims_final.jld2" chlamy_sims
@load "results/chlamy_helices_final.jld2" chlamy_helices
@load "results/phototaxis_pos.jld2"
@load "results/phototaxis_neg.jld2"

f1_beat_plane_orientations = range(0., π/4, 16)
f2_beat_plane_orientations = range(-π/4, π/4, 31)

## Convergence checking
function displacement(N_body=313, Q_body=1213, N_f=27, Q_f=111)
    chlamy_sim = ChlamySimulation("test", beat_plane_tilt=(0.2, 0.0), N_body=N_body, Q_body=Q_body, N_f=N_f, Q_f=Q_f)
    chlamy = generate_chlamy(chlamy_sim)
    prob = SwimmingTrajectoryProblem(chlamy, t_final=1.0, saveat=1.0, eps=0.15)
    solve_problem!(prob)
    prob.traj.x[end]
end

## Visualisation
function animate_sim(sim::ChlamySimulation, N=50; filename=nothing)
    traj50 = continue_periodic_trajectory(sim.traj, N)
    chlamy = generate_chlamy(sim)
    animate(traj50, chlamy, filename=filename, framerate=30)
end


# Paper figures

function plot_chlamy!(parent, sim; rasterize_body=10, kwargs...)
    chlamy = generate_chlamy(sim)

    ax = Axis3(parent[1,1], 
        aspect=:data,
        protrusions=(0,0,0,0)
    )
    hidedecorations!(ax)
    hidespines!(ax)
    viz!(ax, chlamy; rasterize_body=rasterize_body, kwargs...)
    ax
end 

plot_chlamy(sim; kwargs...) = begin
    fig = Figure()
    plot_chlamy!(fig[1,1], sim; kwargs...)
    fig
end


function plot_beat_pattern!(parent, chlamy; 
    num_t=15, 
    show_body=true, 
    bodycolor=Makie.wong_colors()[1], 
    flagella_grad=:greens, 
    kwargs...
)
    g = parent isa GridLayout ? parent : GridLayout(parent)
    ax = Axis3(g[1,1], 
        # title=title,
        aspect=:data, 
        xlabel=L"x \; (\mu\mathrm{m})", 
        ylabel=L"y \; (\mu\mathrm{m})",
        azimuth=-π/2,
        elevation=π/2,
        protrusions=(0,0,0,0)
    )
    # Label(g[1,1,Top()], title)
    hidedecorations!(ax)
    hidespines!(ax)
    location, orientation = chlamy.points.location, chlamy.points.orientation
   
    if show_body
        viz!(ax, chlamy.body, location, orientation; rasterize=10, color=bodycolor)
    end

    cg = cgrad(flagella_grad, num_t)
    ts = range(0.5, 1.5, num_t+1)[1:end-1]

    for (t, color) in zip(ts, cg)
        update_boundary!(chlamy, t)
        f1 = chlamy.flagella[1]
        f2 = chlamy.flagella[2]
        viz!(ax, f1, location, orientation; color=color, kwargs...)
        viz!(ax, f2, location, orientation; color=color, kwargs...)
    end
    ax
end

plot_beat_pattern(chlamy; kwargs...) = begin
    fig = Figure()
    plot_beat_pattern!(fig[1,1], chlamy; kwargs...)
    fig
end

function plot_trajectory!(ax, traj, chlamy; 
    traj_color=:red, 
    traj_linewidth=1.5,
    show_body=true, 
    kwargs...
)
    l = lines!(ax, traj.x, color=traj_color, linewidth=traj_linewidth, colormap=:RdBu, colorrange=(0.0, π))
    if show_body
        move_boundary!(chlamy, traj, length(traj.t))
        viz!(ax, chlamy; rasterize_body=10, kwargs...)
    end
    l   
end

function plot_lambda1_trajectories!(parent, trajs, chlamys;
    show_body=true,
    limits=(-1.0, 85., -20., 60., -10., 60.),
    azimuth=6.952196969386569, 
    elevation=0.3660323641694269,
    colors=Makie.wong_colors()[4:7],
    kwargs...
)
    g = GridLayout(parent[1,1])
    colsize!(g, 1, Relative(1))
    ax = Axis3(g[1, 1], 
        aspect=:data, 
        # limits=limits, 
        xlabeloffset = 20,
        ylabeloffset = 20,
        zlabeloffset = 20,
        azimuth=azimuth,
        elevation=elevation,
        xlabel=L"x \, (\mu\mathrm{m})",
        ylabel=L"y \, (\mu\mathrm{m})",
        zlabel=L"z \, (\mu\mathrm{m})",
        # left, right, bottom, top
        protrusions=(30, 0, 30, 0)
    )
    hidespines!(ax)
    lines = []
    for i in eachindex(trajs)
        push!(lines, plot_trajectory!(ax, trajs[i], chlamys[i]; traj_color=colors[i], traj_linewidth=1.5, show_body=show_body, kwargs...))
    end
    # axislegend(ax, lines, [L"0.0", L"0.1", L"0.26", L"0.47"], L"\Delta\lambda_2", position=:ct, orientation=:horizontal)
    Legend(g[1,2], lines, [L"0.0", L"0.1", L"0.26", L"0.47"], L"\Delta\lambda_2")
    # Box(g[1,2], color=(:red, 0.5))
    colsize!(g, 2, Auto())
    colgap!(g, 0)
    ax
end

plot_lambda1_trajectories(sims=[16, 18, 21, 25], N=50; kwargs...) = begin
    chlamys = [generate_chlamy(chlamy_sims[1,i]) for i in sims]
    trajs = [continue_periodic_trajectory(chlamy_sims[1,i].traj, N) for i in sims]
    fig = Figure()
    plot_lambda1_trajectories!(fig[1,1], trajs, chlamys; kwargs...)
    fig
end


function plot_helix_measurements!(parent)
    g = parent isa GridLayout ? parent : GridLayout(parent)

    lambda2s = f2_beat_plane_orientations
    lambda1s = f1_beat_plane_orientations

    lambda1_labels = [L"0", L"\frac{\pi}{8}", L"\frac{\pi}{4}"]
    lambda2_labels = [L"-\frac{\pi}{4}", L"-\frac{\pi}{8}", L"0", L"\frac{\pi}{8}", L"\frac{\pi}{4}"]

    measurements = [
        axis_velocity,
        curvature,
        axis_polar_angle,
        axis_angular_velocity,
        torsion,
        axis_azimuthal_angle
    ]

    titles = [
        "velocity",
        "curvature",
        "polar direction",
        "angular velocity",
        "torsion",
        "azimuthal direction"
    ]


    labels = [
        L"v\;(\mu\mathrm{m}/\text{beat})",
        L"\kappa\;(\mu\mathrm{m}^{-1})",
        L"\theta\;\text{(rad)}",
        L"\omega\;\text{(rad/beat)}",
        L"\tau\;(\mu\mathrm{m}^{-1})",
        L"\phi\;\text{(rad)}"
    ]

    colormaps = [
        :ice,
        :ice,
        :RdBu,
        Reverse(:cork),
        Reverse(:cork),
        Reverse(:cork)
    ]

    colorranges = [
        (0., 1.25),
        (0., 0.5),
        (π/2-0.4, π/2+0.4),
        (-0.3, 0.3),
        (-0.2, 0.2),
        (-π/2, π/2)
    ]

    colorticks = [
        0:0.5:1,
        0:0.25:0.5,
        ([3π/8, π/2, 5π/8], [L"\frac{3\pi}{8}",L"\frac{\pi}{2}", L"\frac{5\pi}{8}"]),
        -0.2:0.2:0.2,
        -0.2:0.2:0.2,
        (-π/2:π/2:π/2, [L"-\frac{\pi}{2}",L"0", L"\frac{\pi}{2}"])
        ]
    k = 1

    for i in 1:2
        for j in 1:3

            ax = Axis(g[i,2j-1],
                xlabel = i == 2 ? L"\Delta\lambda_2" : "",
                ylabel = j == 1 ? L"\Delta\lambda_1" : "",
                xticks = (-π/4:π/8:π/4, lambda2_labels),
                yticks = (0:π/8:π/4, lambda1_labels),
                yticklabelsvisible = j == 1,
                xticklabelsvisible = i == 2,
                title = titles[k],
                # aspect = 2
            )
            # colsize!(g, 2j-1, Aspect(1.0, 2))

            data = measurements[k].(chlamy_helices)'

            kwargs = colorranges[k] === nothing ?
                (; colormap = colormaps[k]) :
                (; colormap = colormaps[k], colorrange = colorranges[k])

            hm = heatmap!(ax, lambda2s, lambda1s, data; rasterize=10, kwargs...)

            Colorbar(g[i,2j], hm, width = 5, ticks=colorticks[k])
            Label(g[i,2j, Top()], labels[k])

            k += 1
        end
    end

    colgap!(g, 1, 5)
    colgap!(g, 2, 8)
    colgap!(g, 3, 5)
    colgap!(g, 4, 8)
    colgap!(g, 5, 5)
    rowgap!(g, 10)
end

plot_helix_measurements() = begin
    fig = Figure()
    plot_helix_measurements!(fig[1,1])
    fig
end


function beat_plane_phase_portrait!(parent; kwargs...)
    g = parent isa GridLayout ? parent : GridLayout(parent)

    xmin, xmax = -π/4, π/4
    ymin, ymax = -π/4, π/4

    ax = Axis(g[1,1],
        xticks = ([xmin, 0, xmax], [L"-\frac{\pi}{4}", L"0", L"\frac{\pi}{4}"]),
        yticks = ([ymin, 0, ymax], [L"-\frac{\pi}{4}", L"0", L"\frac{\pi}{4}"]),
        aspect = 1,
        xlabel = L"\Delta\lambda_2",
        ylabel = L"\Delta\lambda_1",
        limits = (xmin, xmax, ymin, ymax)
    )
    colsize!(g, 1, Aspect(1.0, 1))

    # shaded regions separated by y = x
    cw_color  = (:yellow, 0.2) # (RGBf(0.82, 0.84, 0.88), 0.35)   # cool slate
    ccw_color = (:orchid2, 0.2) # RGBf(0.89, 0.84, 0.80)   # warm stone

    hidespines!(ax)
    hlines!(ax, [0.0], color = :black, linewidth=0.5)
    vlines!(ax, [0.0], color = :black, linewidth=0.5)
    hidedecorations!(ax, ticklabels = false, label = false)

    lines!(ax, [Point2f(xmin, ymin), Point2f(xmax, ymax)],
        linestyle = :dash, linewidth = 0.5, color = :black)

    lines!(ax, [Point2f(xmax, ymin), Point2f(xmin, ymax)],
        linestyle = :dash, linewidth = 0.5, color = :black)

    poly!(ax, Point2f[
        (xmin, ymin),
        (xmax, ymin),
        (xmax, ymax)
    ], color = cw_color, strokewidth = 0)
    
    poly!(ax, Point2f[
        (xmin, ymin),
        (xmin, ymax),
        (xmax, ymax)
    ], color = ccw_color, strokewidth = 0)
        

    cs = Makie.wong_colors()[4:7]
    lambda1s = [0.0, 0.1, 0.26, 0.47]
    for i in eachindex(lambda1s)
        scatter!(ax, Point2f(lambda1s[i], 0.0 ), color=cs[i], markersize=8, strokewidth=1, strokecolor=:black)
    end

    text!(ax, Point2f(pi/8, pi/8), text=L"\text{circular}", align=(:center, :baseline), rotation=pi/4)
    text!(ax, Point2f(3pi/16, 3pi/16), text=L"\text{decreasing }r", align=(:center, :baseline), rotation=pi/4)
    text!(ax, Point2f(-pi/8, -pi/8), text=L"\text{circular}", align=(:center, :baseline), rotation=pi/4)
    text!(ax, Point2f(-3pi/16, -3pi/16), text=L"\text{decreasing }r", align=(:center, :baseline), rotation=pi/4)
    text!(ax, Point2f(-pi/8, 3pi/16), text=L"\text{helical - increasing }\tau", align=(:left, :center))
    text!(ax, Point2f(-pi/8, pi/8), text=L"\text{linear}", align=(:center, :baseline), rotation=-pi/4)
    text!(ax, Point2f(-3pi/16, 3pi/16), text=L"\text{increasing }\omega", align=(:center, :baseline), rotation=-pi/4)
    text!(ax, Point2f(pi/8, pi/16), text=L"\text{CW helices}", align=(:right, :center))
    text!(ax, Point2f(pi/8, -pi/16), text=L"\text{CCW helices}", align=(:right, :center))
    text!(ax, Point2f(pi/8, -pi/8), text=L"\text{linear}", align=(:center, :baseline), rotation=-pi/4)
    text!(ax, Point2f(3pi/16, -3pi/16), text=L"\text{increasing }\omega", align=(:center, :baseline), rotation=-pi/4)

    ax
end

beat_plane_phase_portrait() = begin
    fig = Figure()
    beat_plane_phase_portrait!(fig[1,1])
    fig
end

function phototaxis_plot!(parent; bodycolor=Makie.wong_colors()[1], color=:greens, show_body=true, kwargs...)
    g = parent isa GridLayout ? parent : GridLayout(parent)

    ax = Axis3(g[1:3,1], 
        aspect=:data,
        # left, right, bottom, top
        protrusions=(30,10,10,0),
        xlabeloffset=20,
        ylabeloffset=20,
        zlabeloffset=20,
        xlabel=L"x \, (\mu\mathrm{m})",
        ylabel=L"y \, (\mu\mathrm{m})",
        zlabel=L"z \, (\mu\mathrm{m})"
    )
    hidespines!(ax)
    cs = Makie.wong_colors()[4:6]

    probs = [prob_rev2, prob_rev]
    chlamy = generate_chlamy(chlamy_sims[1,16])
    modulate_orientation(chlamy, chlamy_sims[1,16].polar_orientation_angle, -0.3, 0.)
    
    mod = [bp_mod_rev2, bp_mod_rev]
    ls = []
    ax2 = Axis(g[3,2], 
        xlabel = L"t"
    )
    hidespines!(ax2)
    
    for i in 1:2
        theta = acos.(getindex.(probs[i].traj.x, 3) ./ norm.(probs[i].traj.x))
        plot_trajectory!(ax, probs[i].traj, chlamy; 
            traj_color=theta, 
            show_body=show_body, 
            traj_linewidth=1.5,
            kwargs...
        )

        βmod = getindex.(mod[i].saveval, 3)
        push!(ls, lines!(ax2, mod[i].t, βmod, color=theta, colormap=:RdBu, colorrange=(0.0, π)))
    end
    Label(g[1:3,1], L"\mathbf{I}", 
        tellwidth=false, tellheight=false,
        halign=:center, valign=:top
    )

    cs = Makie.wong_colors()[4:5]
    move_boundary!(chlamy, [0., 0., 0.], rotation_matrix([0., 0., 1.], -pi/4)*rotation_matrix([1., 0., 0.], 0.2)*rotation_matrix([0., 1., 0.], pi/2-0.3)*rotation_matrix([0., 0., 1], pi/2), 0.0)
    ax3 = plot_beat_pattern!(g[1:2,2], chlamy;
        bodycolor=bodycolor, 
        show_body=show_body, 
        flagella_grad=color, 
        linewidth=2
    )

    # beat planes
    S0 = [Point3(0), Point3(0, -5, 0), Point3(5,-5,0), Point3(5,0,0)]

    B = chlamy.points.orientation
    L1, O1 = chlamy.flagella[1].points.location, chlamy.flagella[1].points.orientation
    L2, O2 = chlamy.flagella[2].points.location, chlamy.flagella[2].points.orientation
    
    S1 = [B*(L1 + O1*p) for p in S0]
    S2 = [B*(L2 + O2*p) for p in S0]
    poly!(ax3, S1, rasterize=10)
    poly!(ax3, S2, rasterize=10)

    arrows3d!(ax3, Point3(B*L1), Point3(B*[1., 0., 0.]), lengthscale=3, rasterize=10)
    arrows3d!(ax3, Point3(B*L2), Point3(B*[1., 0., 0.]), lengthscale=3, rasterize=10)
    
    rowgap!(g, 10)
end 

phototaxis_plot(;show_body=true, kwargs...) = begin
    fig = Figure()
    phototaxis_plot!(fig[1,1], show_body=show_body)
    fig
end