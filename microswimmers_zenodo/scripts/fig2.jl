include("theme.jl")
include("chlamy_analysis.jl")
CairoMakie.activate!()


# a 
chlamy1 = generate_chlamy(chlamy_sims[1,16])
update_boundary!(chlamy1, 0.5)

# b
chlamy2 = generate_chlamy(chlamy_sims[4,25])
update_boundary!(chlamy2, 0.5)
move_boundary!(chlamy2, [0., 0., 0.], rotation_matrix([1., 0., 0.], -π/2)*rotation_matrix([0., 0., 1.], -π/2), 0.0)

# c
sims = [16, 18, 21, 25]
lambda1_chlamys = [generate_chlamy(chlamy_sims[1,i]) for i in sims]
lambda1_trajs = [continue_periodic_trajectory(chlamy_sims[1,i].traj, 50) for i in sims]


function chlamy_layout(;draft=false)
    fig = Figure(size = (APS_FULLWIDTH_PT, 540), padding=(0., 0., 0., 0.))
    g = fig[1,1] = GridLayout()
    
    ga = g[1,1] = GridLayout()
    gb = g[2,1] = GridLayout()
    gc = g[3,1] = GridLayout()
    
    # modelling
    rowsize!(g, 1, Relative(0.3))
    rowsize!(g, 2, Relative(0.25))

    show_body = draft ? false : true
    # x-y view
    ax11 = plot_beat_pattern!(ga[1,1], chlamy1;
        bodycolor=bc, 
        show_body=show_body, 
        flagella_grad=fg, 
        linewidth=2
    )
    
    # scale bar
    lines!(ax11, [Point3(-5.0, -12., 0.0), Point3(0., -12., 0.)], linewidth=3, color=:black)
    
    # text annotations
    t1 = [L"\phi", L"5\mu\mathrm{m}", L"x", L"y", L"z"]
    locs = [Point3(13.0, y, 0.) for y in range(-10., 10., 5)]
    text!(ax11, locs, text=t1, align=(:right, :center))
    # y-z view
    ax12 = plot_beat_pattern!(ga[1,2], chlamy2;
        bodycolor=bc, 
        show_body=show_body, 
        flagella_grad=fg, 
        linewidth=2
    )

     # scale bar
    lines!(ax12, [Point3(-5.0, -6., 0.0), Point3(0., -6., 0.)], linewidth=3, color=:black)
    
      # text annotations
    t2 = [L"\Delta\lambda_1", L"\Delta\lambda_2", L"\theta", L"y", L"z"]
    text!(ax12, locs, text=t2, align=(:right, :center))

    # trajectory plot
    ax13 = plot_lambda1_trajectories!(ga[1,3], lambda1_trajs, lambda1_chlamys;
        bodycolor=bc, 
        show_body=show_body, 
        color=fc, 
        linewidth=0.5
    )
    text!(ax13, Point3(60, 0., 40.), text=L"\Delta\lambda_1=0")

    colsize!(ga, 1, Relative(0.25))
    colsize!(ga, 2, Relative(0.3))
    colgap!(ga, 10)
    
    # heatmaps
    plot_helix_measurements!(gb)
    
    # phase portrait
    beat_plane_phase_portrait!(gc[1,1]; 
        bodycolor=bc, 
        show_body=show_body, 
        color=fg, 
        linewidth=0.5
    )
    colsize!(gc, 1, Relative(0.5))
    
    phototaxis_plot!(gc[1,2]; bodycolor=bc, show_body=show_body, color=fg, linewidth=0.4)
   
    colgap!(gc, 15)
    
    rowgap!(g, 1, 5)
    rowgap!(g, 2, 5)
    fig
end

# generates the trajectories on the beat plane phase portrait
function phase_portrait_examples(;show_body=true)
    fig = Figure(size = (APS_FULLWIDTH_PX // 2, 720))

    centred_trajs = [centred_trajectory(traj) for traj in trajs_pp]

    offset = 10.0
    xs = reduce(vcat, [traj.x for traj in centred_trajs])
    xmin, xmax = extrema(p[1] for p in xs)
    ymin, ymax = extrema(p[2] for p in xs)
    zmin, zmax = extrema(p[3] for p in xs)
    limits=(xmin - offset, xmax + offset, ymin - offset, ymax + offset, zmin - offset, zmax + 5.0)

    azimuth = 5.5
    elevation = 0 # 0.0926990816987243
    for i in eachindex(chlamys_pp)

        ax = Axis3(fig[i,1],
            aspect = :data,
            protrusions = (0,0,0,0),
            # viewmode = :fit,
            backgroundcolor = :transparent,
            azimuth=azimuth,
            elevation=elevation,
            limits=limits
        )
    
        hidedecorations!(ax)
        hidespines!(ax)
    
        plot_trajectory!(ax, centred_trajs[i], chlamys_pp[i];
            traj_linewidth = 1.0,
            show_body = show_body,
            bodycolor=bc,
            color=fc
        )
    end
    fig
end


# make the figure
fig = chlamy_layout(draft=false)