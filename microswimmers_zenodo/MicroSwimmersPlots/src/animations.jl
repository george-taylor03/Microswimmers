function animate(microswimmer::MicroSwimmer, T=5.0, num_t=20*5+1)
    fig = Figure()
    ax = Axis3(fig[1,1], aspect=:data)
    B = viz!(ax, microswimmer)
    display(fig)
    for t in range(0, T, num_t)
        update_boundary!(microswimmer, t)
        update_buffer_observable!(B, microswimmer)
        sleep(0.01)
    end
end

function animate(
    traj::Trajectory,
    microswimmer::MicroSwimmer;
    wall=false, 
    limits=nothing, 
    step=5,
    filename=nothing,
    framerate=30,
    compression=20,
    elevation=π/36, 
    azimuth=π/4,
    azimuth_rate=nothing,
    fig_size=(1920, 1080),
    ax_visible=true,
    n_fade=1000,
    kwargs...
)
    ts = traj.t
    
    t = Observable(1)
    T = @lift(traj.x[1:$t])

    trail_color = @lift begin
        n = length($T)
        start = max(1, n - n_fade)
        # alpha ramps 0→1 over the tail
        [RGBAf(1, 0, 0, (i - start) / max(1, n - start)) for i in start:n]
    end
    trail_pts = @lift($T[max(1, length($T) - n_fade):end])

    
    if isnothing(limits)
        r_max = maximum(norm, microswimmer.points.force_pts)
        xs = getindex.(traj.x, 1)
        ys = getindex.(traj.x, 2)
        zs = getindex.(traj.x, 3)
        limits = (
            (minimum(xs) - r_max, maximum(xs) + r_max),
            (minimum(ys) - r_max, maximum(ys) + r_max),
            (minimum(zs) - r_max, maximum(zs) + r_max)
        )
    end
        
    fig = Figure(size=fig_size)
    ax = Axis3(fig[1, 1], aspect=:data, limits=limits, elevation=elevation, azimuth=azimuth)
    if !ax_visible
        hidedecorations!(ax)
        hidespines!(ax)
    end
    obs = viz!(ax, microswimmer; kwargs...)
    lines!(ax, trail_pts, color=trail_color, linewidth=4)
    # lines!(ax, T, color=:red, linewidth=0.5)
    yield()
    display(fig)
    
    on(t) do i
        if !isnothing(azimuth_rate)
            ax.azimuth[] = azimuth + i * azimuth_rate
        end
        move_boundary!(microswimmer, traj.x[i], traj.b1[i], traj.b2[i], ts[i])
        update_buffer_observable!(obs, microswimmer)
    end

    if wall
        v = [
            -1. -1. 0.;
            -1.  1. 0.;
             1.  1. 0.;
             1. -1. 0.
        ]
        f = [1 2 3; 3 4 1]
        mesh!(ax, v, f)
    end

    if isnothing(filename)
        for i in 1:step:length(ts)
            t[] = i
            # yield()
            sleep(1 // framerate)
        end
    else
        record(fig, filename, 1:step:length(ts); framerate=framerate, compression=compression) do i 
            t[] = i
        end 
    end
    # fig
end

animate(prob::SwimmingTrajectoryProblem; kwargs...) = animate(
        prob.traj, prob.swimming_problem.microswimmer; 
        kwargs...   
)

function animate!(
    ax, t::Observable{Int},
    traj::Trajectory, microswimmer::MicroSwimmer;
    n_fade=60, trail_color=:red, kwargs...
)

    lo = fill(Inf, 3); hi = fill(-Inf, 3)
    for i in eachindex(traj.x)
        move_boundary!(microswimmer, traj.x[i], traj.b1[i], traj.b2[i], traj.t[i])
        for p in eachcol(microswimmer.points.force_pts)   # now in world frame after move_boundary!
            lo .= min.(lo, p .+ traj.x[i])
            hi .= max.(hi, p .+ traj.x[i])
        end
    end
    pad = 0.05 .* (hi .- lo)   # small margin
    ax.limits = (
        (lo[1]-pad[1], hi[1]+pad[1]),
        (lo[2]-pad[2], hi[2]+pad[2]),
        (lo[3]-pad[3], hi[3]+pad[3]),
    )
    @info "limits" ax.limits[]


    ts = traj.t
    T = @lift(traj.x[1:$t])

    obs = viz!(ax, microswimmer; kwargs...)

    # two-layer fading trail
    trail_pts = @lift($T[max(1, length($T) - n_fade):end])
    lines!(ax, T, color=(trail_color, 0.15), linewidth=1)
    lines!(ax, trail_pts, color=trail_color, linewidth=4)

    on(t) do i
        move_boundary!(microswimmer, traj.x[i], traj.b1[i], traj.b2[i], ts[i])
        update_buffer_observable!(obs, microswimmer)
    end
    ax
end

function animate_panels(
    trajs::Vector{<:Trajectory}, swimmers::Vector{<:MicroSwimmer};
    titles=nothing, filename=nothing, framerate=30, compression=20,
    elevation=π/36, azimuth=π/4, fig_size=(1920, 720), step=2,
    limits=nothing, kwargs...
)
    n = length(trajs)
    fig = Figure(size=fig_size)
    display(fig)
    t = Observable(1)
    for k in 1:n
        ax = Axis3(fig[1, k], aspect=:data,
                #    limits = isnothing(limits) ? fill(nothing, 6) : limits,
                   elevation=elevation, azimuth=azimuth)
        isnothing(titles) || (ax.title = titles[k])
        animate!(ax, t, trajs[k], swimmers[k]; kwargs...)
    end

    n_frames = maximum(length(tr.t) for tr in trajs)

    if isnothing(filename)
        for i in 1:step:n_frames
            t[] = i
            sleep(1 // framerate)
        end
    else
        record(fig, filename, 1:step:n_frames; framerate=framerate, compression=compression) do i
            t[] = i
        end
    end
end


"""Animate particle trajectories in a flow, NEEDS UPDATING"""
function animate(
    ts::Vector,
    traj::Matrix,
    microswimmer::MicroSwimmer;
    wall=false, 
    limits=nothing, 
    step=5,
    traj_length=600,
    filename=nothing,
    framerate=30,
    elevation=π/6, 
    azimuth=π/4,
    num_particles=nothing,
    kwargs...
)
    t = Observable(1)
    num_particles = isnothing(num_particles) ? size(traj,1) ÷ 3 : num_particles

    particle_trajectories = [@lift(traj[3i-2:3i, max(1, $t-traj_length):$t]) for i in 1:num_particles]

    if isnothing(limits)
        r_max = maximum(norm, microswimmer.points.force_pts)
        limits = (
            (- 1.4r_max, 1.4r_max),
            (- 1.4r_max, 1.4r_max),
            (- 1.4r_max, 1.4r_max)
        )
    end

    fig = Figure()
    ax = Axis3(fig[1, 1], aspect=:data, limits=limits, elevation=elevation, azimuth=azimuth)
    obs = viz!(ax, microswimmer; kwargs...)


    for i in 1:num_particles
        lines!(ax, particle_trajectories[i], color=RGBf(0.2,0.6,0.6), linewidth=0.7)
        scatter!(ax, @lift(Point3($(particle_trajectories[i])[:,end])), color=RGBf(0.2,0.6,0.6), markersize=5)
    end

    if wall
        v = [
            -1. -1. 0.;
            -1.  1. 0.;
                1.  1. 0.;
                1. -1. 0.
        ]
        faces = [1 2 3; 3 4 1]
        mesh!(ax, v, faces)
    end
    display(fig)

    
    on(t) do i
        update_boundary!(microswimmer, ts[i])
        update_buffer_observable!(obs, microswimmer)
    end
    if isnothing(filename)
        for i in 1:step:length(ts)
            t[] = i
            sleep(0.05)
        end
    else
        record(fig, filename, 1:length(ts); framerate=framerate) do i 
            t[] = i
        end 
    end
end

animate(prob::ParticleTrajectoryProblem; kwargs...) = animate(
        prob.t, prob.trajectories, prob.resistance_problem.microswimmer; 
        kwargs...   
)