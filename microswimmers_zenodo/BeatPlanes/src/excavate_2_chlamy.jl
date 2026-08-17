mutable struct Excavate2ChlamySimulation <: Simulation
    # Metadata
    sim_type::String             # e.g., "type1_centerline"
    body_dims::Tuple{Float64, Float64, Float64}
    model_posterior::PlanarFlagellum
    model_anterior::PlanarFlagellum

    attachment_polar_angle::Float64       # considered the same for both flagella
    attachment_azimuthal_angle::Float64  
    groove_height::Float64
    posterior_angle_about_y::Float64
    anterior_angle_about_x::Float64
    posterior_curvature::Float64

    N_body::Int
    Q_body::Int
    N_f::Int
    Q_f::Int
    eps::Float64

    # Outputs
    traj::Union{Nothing, Trajectory{Float64}}
end

function Excavate2ChlamySimulation(
    type::String;
    a=4.0,
    b=3.0,
    c=3.0,
    L=10.,
    α=π/7,   # attachment angle of flagella in x-y plane relative to x-axis (polar angle)
    β=π/2,    # attachment angle of flagella in y-z plane relative to y-axis (azimuthal angle)
    d=2.,
    Ry_1=-8π/9,
    Rx_2=-π/2,
    C=0.,
    N_body=234,
    Q_body=1826,
    N_f=54,
    Q_f=234,
    eps=0.1
)
    Excavate2ChlamySimulation(
        type,
        (a,b,c),
        PlanarFlagellum(L, C, 0.7, 0.15, 2π, 2π, 2π, 0.0),
        PlanarFlagellum(L, -2.5, 0.7, 0.15, 2π, 2π, 2π, 0.0),
        α, β,
        d,
        Ry_1, Rx_2,
        C,
        N_body, Q_body, N_f, Q_f,
        eps,
        nothing
    )
end

function generate_e2c(sim::Excavate2ChlamySimulation)
    a, b, c = sim.body_dims
    d = sim.groove_height
    α = sim.attachment_polar_angle
    β = sim.attachment_azimuthal_angle
    Ry_1 = sim.posterior_angle_about_y
    Rx_2 = sim.anterior_angle_about_x

    body = CellBody(
        EllipsoidalGroovedBody(a, b, c, [0., 0., d]),
        sim.N_body, sim.Q_body
    )
    
    posterior = Flagellum(
        sim.model_posterior, sim.N_f, sim.Q_f,
        location=[a*cos(α), b*sin(α)*cos(β), b*sin(α)*sin(β)],
        orientation=rotation_matrix([0.,1.,0.], Ry_1)*rotation_matrix([1.,0.,0.], 1.0π)
    )
    anterior = Flagellum(
        sim.model_anterior, sim.N_f, sim.Q_f,
        location=[a*cos(α), -b*sin(α)*cos(β), -b*sin(α)*sin(β)]  ,
        orientation=rotation_matrix([1.,0.,0.], Rx_2)
    )
    
    Flagellate(body, [posterior, anterior])
end

function run!(sim::Excavate2ChlamySimulation)
    e2c = generate_e2c(sim)
    prob = SwimmingTrajectoryProblem(e2c, t_final=1.0, saveat=0.05, eps=sim.eps)
    solve_problem!(prob, periodic=true)
    sim.traj = prob.traj   
end

