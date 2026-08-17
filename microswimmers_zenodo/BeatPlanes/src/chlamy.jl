mutable struct ChlamySimulation{M <: FlagellumModel} <: Simulation
    # Metadata
    sim_type::String           
    body_dims::Tuple{Float64, Float64, Float64}
    model::M
    beat_plane_tilt::Tuple{Float64, Float64}
    polar_location_angle::Float64    # attachment angle of flagella
    polar_orientation_angle::Float64 # tangent angle of flagella at the base

    N_body::Int
    Q_body::Int
    N_f::Int
    Q_f::Int
    eps::Float64

    # Outputs 
    traj::Union{Nothing, Trajectory{Float64}}
end

function ChlamySimulation(
    type::String;
    # StandingWaveFlagellum parameters
    L=11.0,
    C=-2.2,
    A01=0.5976685766851538,
    ϕ01=-0.9891230066424034,
    A11=0.5367767361004127,
    ϕ11=0.5651694693686302,
    A21=0.15781298240563,
    ϕ21=2.1971231087346057,
    A31=0.0,
    ϕ31=0.0,
    ω=2π,
    # Simulation parameters
    beat_plane_tilt=(0.0, 0.0),
    α=0.3,                  # polar_location_angle
    γ=0.5,                  # polar_orientation_angle
    a=5.0, b=4.0, c=4.0,
    N_body=237, Q_body=999,
    N_f=87, Q_f=261,
    eps=0.15
)
    model = StandingWaveFlagellum(L, C, A01, ϕ01, A11, ϕ11, A21, ϕ21, A31, ϕ31, ω)
    ChlamySimulation(
        type,
        (a, b, c),
        model,
        beat_plane_tilt,
        α, γ,
        N_body, Q_body, N_f, Q_f,
        eps,
        nothing
    )
end

# Use for an arbitrary flagellum model
function ChlamySimulation(
    type::String,
    model::FlagellumModel;
    α=0.3, 
    Δλ₁=0.0, Δλ₂=0.0,
    γ=0.5,
    a=5.0, b=4.0, c=4.0,
    N_body=713, Q_body=4317,
    N_f=27, Q_f=117,
    eps=0.1
)
    ChlamySimulation(
        type,
        (a,b,c),
        model,
        (Δλ₁, Δλ₂),
        α, γ,
        N_body, Q_body, N_f, Q_f,
        eps,
        nothing
    )
end

function change_disc(sim::ChlamySimulation, N_body, Q_body, N_f, Q_f)
    sim.N_body = N_body
    sim.Q_body = Q_body
    sim.N_f = N_f
    sim.Q_f = Q_f
end

function generate_chlamy(sim::ChlamySimulation)
    a,b,c = sim.body_dims
    model = sim.model
    α = sim.polar_location_angle
    Δλ₁, Δλ₂ = sim.beat_plane_tilt
    γ = sim.polar_orientation_angle


    body = CellBody(
        EllipsoidBody(a, b, c), 
        sim.N_body, sim.Q_body
    )
    f1 = Flagellum(
        model, sim.N_f, sim.Q_f,     
        location=[a*cos(α), b*sin(α), 0.0],
        orientation=rotation_matrix([0.0, 0.0, 1.0], γ)*rotation_matrix([1.0, 0.0, 0.0], π-Δλ₁) 
    )

    f2 = Flagellum(
        model, sim.N_f, sim.Q_f,
        location=[a*cos(α), -b*sin(α), 0.0],
        orientation=rotation_matrix([0.0, 0.0, 1.0], -γ)*rotation_matrix([1.0, 0.0, 0.0], Δλ₂),
    )
    Flagellate(body, [f1, f2])
end

function run!(sim::ChlamySimulation)
    chlamy = generate_chlamy(sim)
    prob = SwimmingTrajectoryProblem(chlamy, t_final=1.0, saveat=0.05, eps=sim.eps)
    solve_problem!(prob, periodic=true)
    sim.traj = prob.traj
end


# Phototaxis
eyespot_orientation = @SVector [0.0, 1.0/sqrt(2.0), -1.0/sqrt(2.0)]
ex = @SVector [1.0, 0.0, 0.0]
ey = @SVector [0.0, 1.0, 0.0]
ez = @SVector [0.0, 0.0, 1.0]
heaviside(x) = x ≥ 0 ? 1.0 : 0.0

function beatplane_modulation(I0, σ, B)
    N = B*eyespot_orientation
    alignment = dot(N, ez)
    σ*log(1.0 + I0*max(0.0, alignment))
end

function modulate_orientation(chlamy, γ, β, bp_mod)
    chlamy.flagella[1].points.orientation = rotation_matrix(ez, γ)*rotation_matrix(ex, -bp_mod)*rotation_matrix(ey, -β)*rotation_matrix(ex, π)
    chlamy.flagella[2].points.orientation = rotation_matrix(ez, -γ)*rotation_matrix(ex, bp_mod)*rotation_matrix(ey, β)
end 

function PhototaxisProblem(
    sim::ChlamySimulation;
    x0_0=SVector(0.0, 0.0, 0.0),
    b1_0=SVector(1.0, 0.0, 0.0),
    b2_0=SVector(0.0, 1.0, 0.0),  
    t_final=20.0,
    saveat=0.05,
    eps=0.01,
    mu=1.0,
    I0=20.0,
    σ=1.0,
    β=-0.3
)
    γ = sim.polar_orientation_angle
    S = generate_chlamy(sim)

    T = eltype(S.points.force_pts)
    # x0 = SVector{3,T}(x0)
    # B = SMatrix{3,3,T}(B)

    sprob = SwimmingProblem(S; eps=T(eps), mu=T(mu))
    # es = B*eyespot_orientation

    # x0_0 = SVector{3,T}(0, 0, 0)
    # b1_0 = SVector{3,T}(1, 0, 0)
    # b2_0 = SVector{3,T}(0, 1, 0)
    X0 = vcat(x0_0, b1_0, b2_0)

    function rhs(X, p, t)
        x0 = SVector{3,T}(X[1:3])
        b1 = SVector{3,T}(X[4:6])
        b2 = SVector{3,T}(X[7:9])

        B = hcat(b1, b2, cross(b1,b2))
        bp_mod = beatplane_modulation(I0, σ, B)
        modulate_orientation(S, γ, β, bp_mod)

        move_boundary!(sprob, x0, b1, b2, t)
        solve_problem!(sprob)
        Ω = get_Ω(sprob)
        vcat(get_U(sprob), cross(Ω, b1), cross(Ω, b2))
    end

    function save_beatplane_tilt(X, t, integ)
        b1 = SVector{3,T}(X[4:6])
        b2 = SVector{3,T}(X[7:9])
    
        B = hcat(b1, b2, cross(b1,b2))
        N = B*eyespot_orientation
        alignment = -dot(N, -ez)
        intensity = I0*alignment*heaviside(alignment)
        bp_mod = σ*log(1.0 + intensity)
        # beatplane_modulation(I0, σ, B)
        (alignment, intensity, bp_mod)
    end
    
    beatplane_values = SavedValues(Float64, Tuple{Float64, Float64, Float64})
    cb = SavingCallback(save_beatplane_tilt, beatplane_values, saveat=T(saveat))

    SwimmingTrajectoryProblem(
        sprob,
        ODEProblem(
            rhs, 
            X0, 
            (T(0), T(t_final)), 
            saveat=T(saveat),
            callback=cb
        ),
        nothing
    ), beatplane_values
end 

function eyespot_location(sim::ChlamySimulation)
    a,b,c = sim.body_dims
    @SVector [0.0, b*cos(π/4), -b*sin(π/4)]
end

function eyespot(chlamy::MicroSwimmer)
    loc = chlamy.points.location
    or = chlamy.points.orientation
    b = chlamy.body.model.b
    @info "" or 
    loc + or*[0.0, b*cos(π/4), -b*sin(π/4)], or*eyespot_orientation
end

eyespot(prob::SwimmingTrajectoryProblem) = eyespot(prob.swimming_problem.microswimmer)


