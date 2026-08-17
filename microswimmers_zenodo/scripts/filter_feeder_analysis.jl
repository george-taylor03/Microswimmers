using BeatPlanes
using JLD2
using MicroSwimmers
using MicroSwimmersPlots
using GLMakie


@load "results/filter_feeder_sims_final.jld2" ff_sims
@load "results/filter_feeder_inverted_final.jld2" inverted_sims
@load "results/sphere_flow.jld2"
@load "results/cavity_flow.jld2"


function animate_sim(sim::FilterFeederSimulation; filename=nothing)
    ff = generate_filter_feeder(sim)
    animate(sim.t, sim.traj, ff, limits=(-31., 25, -25, 25, -25, 25), filename=filename, framerate=30)
end


# for convergence checking
function flux_and_power(N_body=313, Q_body=1213, N_f=17, Q_f=111)
    sim = FilterFeederSimulation("test", N_body=N_body, Q_body=Q_body, N_f =N_f, Q_f =Q_f, h=-9.0, waves=20, eps=0.15)
    run!(sim, t_final=0.0)
    sim.flux, sim.power
end


# paper figures
function plot_beat_pattern!(parent, ff; 
    num_t=15, 
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
    location, orientation = ff.points.location, ff.points.orientation

    cg = cgrad(flagella_grad, num_t)
    ts = range(0.5, 1.5, num_t+1)[1:end-1]

    for (t, color) in zip(ts, cg)
        update_boundary!(ff, t)
        f1 = ff.flagella[1]
        viz!(ax, f1, location, orientation; color=color, kwargs...)
    end
    ax
end

function plot_particle_trajectories(;kwargs...)
    fig = Figure()
    plot_particle_trajectories!(fig[1,1], generate_filter_feeder(ff_sim), trajs; kwargs...)
    fig
end

function plot_particle_trajectories!(parent, ff, trajs; kwargs...)
    g = parent isa GridLayout ? parent : GridLayout(parent)
    
    N = length(trajs[1])
    alphas = range(0.1, 0.8, length=N)
    colors = [RGBAf(0.2, 0.6, 0.6, a) for a in alphas]

    ax = Axis3(g[1,1], aspect=:data, protrusions=(0,0,0,0))
    hidedecorations!(ax)
    hidespines!(ax)

    viz!(ax, ff; kwargs...)
    lines!.(ax, trajs, color=RGBf(0.2,0.6,0.6), linewidth=0.5)
end

function plot_filter_feeder!(parent, sim; rasterize_body=1, show_feeding_surface=true, kwargs...)
    ff = generate_filter_feeder(sim)

    ax = Axis3(parent[1,1], 
        aspect=:data,
        protrusions=(0,0,0,0)
    )
    hidedecorations!(ax)
    hidespines!(ax)
    viz!(ax, ff; rasterize_body=rasterize_body, kwargs...)
    if show_feeding_surface
        r = sqrt(sim.body_radius^2 - sim.band_height^2)
        t = range(0, 2π, 100)
        xyz = Point3f.(-30., r*cos.(t), r*sin.(t))
        lines!(ax, xyz, color=(:red, 0.4))
    end

    ax
end 

plot_filter_feeder(sim) = begin
    fig = Figure()
    plot_filter_feeder!(fig[1,1], sim)
    fig
end

function efficiency_plot!(parent)
    @info "" parent
    g = parent isa GridLayout ? parent : GridLayout(parent)

    # --- shared data ---
    h = range(0., 19.0, 16)
    hai = h ./ 20.0
    halfds = hai[1]/2

    # --- top: heatmap ---
    ax_top = Axis(g[1,1],
        title="feeding efficiency per cilium",
        ylabel=L"w",
        limits=(hai[1] - halfds, hai[end] + halfds, nothing, nothing)
    )

    r = [sqrt(sim.body_radius^2 - sim.band_height^2) for sim in ff_sims]
    num_flagella = [floor(2π*ri/sim.flagellum_spacing) for (ri,sim) in zip(r, ff_sims)]
    Q_per_cilium = getproperty.(ff_sims, :flux) ./ num_flagella
    P_per_cilium = getproperty.(ff_sims, :power) ./ num_flagella
    eff_per_cilium = Q_per_cilium.^2 ./ P_per_cilium
    not_right = 2.0*[sim.metachronal_waves for sim in ff_sims] .> num_flagella
    eff_per_cilium[not_right] .= NaN

    h1 = heatmap!(ax_top, hai, 0:30, eff_per_cilium, colormap=:matter, rasterize=10)
    hlines!(ax_top, 5:5:25, color=ColorSchemes.okabe_ito[1:5], linestyle=:dash)

    # --- bottom: line plot ---
    h2 = range(0., 19.0, 16)[3:end]
    hai2 = h2 ./ 20
    sphere_sims = ff_sims[3:end, 6:5:26]
    inv_sims = inverted_sims[2:end,:]

    ax_bot = Axis(g[2,1],
        title="cavity vs. sphere efficiency",
        xlabel=L"h/a",
        ylabel=L"\varepsilon",
        limits=(hai[1], hai[end], nothing, nothing)
    )

    cs = ColorSchemes.okabe_ito

    function get_effs(sims)
        r = [sqrt(sim.body_radius^2 - sim.band_height^2) for sim in sims]
        nf = [floor(2π*ri/sim.flagellum_spacing) for (ri,sim) in zip(r,sims)]
        effs = (getproperty.(sims,:flux)./nf).^2 ./ (getproperty.(sims,:power)./nf)
        not_right = 2.0 .* [sim.metachronal_waves for sim in sims] .> nf
        effs[not_right] .= NaN
        effs
    end

    s_effs = get_effs(sphere_sims)
    inv_effs = get_effs(inv_sims)

    for i in 1:5
        lines!(ax_bot, hai2, inv_effs[:,i], color=cs[i], label=string(5i))
        lines!(ax_bot, hai2, s_effs[:,i], color=cs[i], linestyle=:dash)
        for (col, marker) in [(inv_effs[:,i], :circle), (s_effs[:,i], :diamond)]
            if any(isfinite, col)
                idx = argmax(replace(col, NaN => -Inf))
                scatter!(ax_bot, [hai2[idx]], [col[idx]], color=cs[i], marker=marker, markersize=5)
            end
        end
    end

    # --- shared axes ---
    linkxaxes!(ax_top, ax_bot)
    ax_top.xticklabelsvisible = false
    ax_top.xlabelvisible = false

    # --- colorbar and legend, minimal width ---
    Colorbar(g[1,2], h1, width=5, halign=:left)
    Label(g[1,2,Top()], L"\varepsilon")
    Legend(g[2,2], ax_bot, L"w", rowgap=0, framevisible=false, halign=:left)
    # axislegend(ax_bot, L"w", position=:lt, rowgap=0, framevisible=false)
    colgap!(g, 4)
    rowgap!(g, 4)
    colsize!(g, 2, Relative(0.05))

    (ax_top, ax_bot)
end

efficiency_plot() = begin
    fig = Figure()
    efficiency_plot!(fig[1,1])
    fig
end

function flow_plot_comparison!(parent)
    g = parent isa GridLayout ? parent : GridLayout(parent)

    iff = ff_cavity
    in_ff(x) = is_inside_ellipsoid(x, zeros(3), [20,20,20]) && !is_inside_ellipsoid(x, iff.body.model.groove_center, [20,20,20])
    stream_kwargs = (;
        arrow_size=4.,
        linewidth=.8
    )
    stream_kwargs=(;
        gridsize=(128,128,1), 
        density=0.45, 
        arrow_size=3, 
        linewidth=0.6,
        stepsize=1e-2
    )
    crange = (1e-2, 4.55)
    cmap = :viridis
    bc = :gray50

    hm_kwargs = (;
        colormap=cmap, 
        colorscale=log10, 
        colorrange=crange,
        nan_color=bc,
        rasterize=10
    )

    ax1 = stream!(g[1,1], ave_vf; 
        in_domain = (x,z) -> x^2 + z^2 > 20.0^2,
        show_colorbar=false, 
        stream_kwargs=stream_kwargs,
        hm_kwargs=hm_kwargs
    )
    ax1.xticklabelsvisible = false
    ax1.xlabelvisible = false

    # 
    r = 20.0 # sphere radius
    h = -8.866666666666667 # band height
    z = sqrt(r^2 - h^2)
    lines!(ax1, [Point2f(-30,-z), Point2f(-30,z)], color=:red, linewidth=2)
    
    ax2 = stream!(g[2,1], iave_vf; 
        in_domain = (x,z) -> !in_ff([x, 0., z]),
        show_colorbar=false,
        stream_kwargs=stream_kwargs,
        hm_kwargs=hm_kwargs
    )
    lines!(ax2, [Point2f(-30,-z), Point2f(-30,z)], color=:red, linewidth=2)

    for ax in [ax1, ax2]
        ax.xlabelpadding=1
        ax.ylabelpadding=1
    end
    
    Colorbar(g[1,2], 
        limits=crange, 
        colormap=cmap, 
        scale=log10,
        width=5,
    )
    Colorbar(g[2,2], 
        limits=crange, 
        colormap=cmap, 
        scale=log10,
        width=5,
    )
    Label(g[1,2, Top()], L"v \, (\mu\mathrm{m/beat})")
    Label(g[2,2, Top()], L"v \, (\mu\mathrm{m/beat})")

    colsize!(g, 1, Aspect(1, 8/7.5))
    colgap!(g, 0)
    rowgap!(g, 5)
    g
end

flow_plot_comparison() = begin
    fig = Figure()
    flow_plot_comparison!(fig[1,1])
    fig
end

