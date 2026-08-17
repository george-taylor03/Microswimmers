mutable struct FilterFeederSimulation <: Simulation
    # Metadata
    sim_type::String             # e.g., "type1_centerline"
    body_radius::Float64
    band_height::Float64
    flagellum_spacing::Float64
    metachronal_waves::Int64

    model::PlanarFlagellum
    beat_plane_tilt::Float64

    N_body::Int
    Q_body::Int
    N_f::Int
    Q_f::Int
    eps::Float64

    # Outputs
    t::Union{Nothing, Vector{Float64}}
    traj::Union{Nothing, Matrix{Float64}}
    flux::Union{Nothing, Float64}
    power::Union{Nothing, Float64}
end

function FilterFeederSimulation(
    type::String;
    L=11.0,
    C=-2.46329, 
    R₀=0.35, 
    R₁=1.5, 
    k=0.7, 
    ϕ=6.218, 
    ω=2π, δ=0.0,
    β=0.0,
    a=20.0,
    h=-17.0,
    d=2.0,
    waves=6,
    N_body=113, Q_body=517,
    N_f=15, Q_f=117,
    eps=0.1
)

    FilterFeederSimulation(
        type,
        a, h, d, waves,
        PlanarFlagellum(L, C, R₀, R₁, k, ϕ, ω, 0.0),
        β,
        N_body, Q_body, N_f, Q_f,
        eps,
        nothing,
        nothing,
        nothing,
        nothing
    )
end


# adding 'inverted' anywhere in the sim_type string generates the cavitated body 
function generate_body(sim::FilterFeederSimulation)
    a = sim.body_radius
    h = sim.band_height

   if occursin("inverted", sim.sim_type)
        # Inverted mouth geometry – placeholder example.
        r = sqrt(sim.body_radius^2 - sim.band_height^2)

        body_geom = EllipsoidalGroovedBody(
            a, a, a,          # outer ellipsoid radii
            a, a, a,          # inner groove radii
            [2h, 0.0, 0.0],     # groove offset
            orientation=rotation_matrix([0., 1., 0.], -π/2)
        )
        CellBody(body_geom, sim.N_body, sim.Q_body)
    else
        SphericalBody(a; N=sim.N_body, Q=sim.Q_body)
    end
end
    

function generate_filter_feeder(sim::FilterFeederSimulation)
    body = generate_body(sim)
    
    R_tilt = rotation_matrix([1., 0., 0.], sim.beat_plane_tilt)
    R = rotation_matrix([0., 0., 1.0], π/2)
    r = sqrt(sim.body_radius^2 - sim.band_height^2)
    θ = atan(sim.band_height, r)
    num_flagella = floor(2π*r/sim.flagellum_spacing)
    
    a = sim.body_radius
    h = sim.band_height
    
    @unpack L, C, R₀, R₁, k, δ, ω, ϕ = sim.model 
    # generate flagella positioned around the sphere
    flagella = [Flagellum(
        PlanarFlagellum(L, C, R₀, R₁, k, ϕ, ω, 2π*n*sim.metachronal_waves/num_flagella),
        sim.N_f, 
        sim.Q_f, 
        location=[sim.band_height, r*cos(2π*n/num_flagella), r*sin(2π*n/num_flagella)],
        orientation = rotation_matrix([1., 0., 0.], 2π*n/num_flagella) * rotation_matrix([0., 0., 1.], -θ) * R * R_tilt
    )  for n in 1:num_flagella] 

    Flagellate(
        body,
        flagella
    )
end

function generate_filter_feeder_random_phase(sim::FilterFeederSimulation)
    body = generate_body(sim)
    
    R_tilt = rotation_matrix([1., 0., 0.], sim.beat_plane_tilt)
    R = rotation_matrix([0., 0., 1.0], π/2)
    r = sqrt(sim.body_radius^2 - sim.band_height^2)
    θ = atan(sim.band_height, r)
    num_flagella = floor(2π*r/sim.flagellum_spacing)
    
    @unpack L, C, R₀, R₁, k, δ, ω, ϕ = sim.model 

    phases = 2π .* rand(Int(num_flagella))

    flagella = [Flagellum(
        PlanarFlagellum(L, C, R₀, R₁, k, ϕ, ω, phases[n]),
        sim.N_f, 
        sim.Q_f, 
        location=[sim.band_height, r*cos(2π*n/num_flagella), r*sin(2π*n/num_flagella)],
        orientation = rotation_matrix([1., 0., 0.], 2π*n/num_flagella) * rotation_matrix([0., 0., 1.], -θ) * R * R_tilt
    ) for n in 1:Int(num_flagella)]

    Flagellate(body, flagella)
end

# set t_final > 0 if you also want to run the ParticleTrajectoryProblem
function run!(sim::FilterFeederSimulation; x=-25.0, t_final=0.0, saveat=0.05, num_t=10, ff_fn=generate_filter_feeder)
    ff = ff_fn(sim)
    if t_final > 0.0
        prob = ParticleTrajectoryProblem(
            ff, 
            ys=range(-sim.body_radius, sim.body_radius, 10), 
            zs=range(-sim.body_radius, sim.body_radius, 10), 
            x=x,
            t_final=t_final,
            saveat=saveat,
            eps=0.1
        )
        solve_problem!(prob)
        sim.t = prob.t
        sim.traj = prob.trajectories
    end

    # flux and power
    rprob = ResistanceProblem(ff, eps=sim.eps)
    u = FluidVelocity(rprob)
    
    r = sqrt(sim.body_radius^2 - sim.band_height^2)
    z_top(y) = sqrt(r^2 - y^2)
    z_bot(y) = -sqrt(r^2 - y^2)
    
    fluxes = []
    powers = []
    for t in range(0,1,num_t)[1:end-1]
        update_boundary!(rprob, t)
        solve_problem!(rprob)
        push!(fluxes, velocity_flux_polar(u, x, 0.0, 0.0, r))
        push!(powers, total_power(rprob))
    end     
    @info "" fluxes
    @info "" powers
    sim.flux = mean(fluxes)
    sim.power = mean(powers)
    fluxes
end





