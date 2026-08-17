#Imports
using MicroSwimmers
using MicroSwimmersPlots
using GLMakie
include("excavate_body_design.jl")

 
##Design Basic Fallugm

omega_beat = 6.283185307179586

##Iterate Through Beat pattern for 
for mult in 1:5
    println("Beat Frequency: $mult")
    println(omega_beat * mult)

    #Posterior flagellum
    f = PlanarStandingWaveFlagellum{Float64}(10.0, omega_beat * mult, 0, [0.15, 0.0, -0.35, 0.0], [-0.3, 0.4, 0.0, -0.3])

    posterior = PlanarVanedFlagellum(f, 0.1, 0.9, 1.0)

    #anterior flagellum
    anterior = ThreeDimensionalFlagellum{Float64}(9.0, 1.0, 0.0, 1.16, 14.0, 0.16, 1.0, 0.8, 0.53, 21.0, -0.16, 0.0, 0.3584073464102069)

    #Design SuperEllipsoid Body
    # jakoba parameters
    el = SuperEllipsoid(3.9, 2.2, 2.2)
    groove = Posed(SuperEllipsoid(3.9, 2.2, 2.2; κx = 0.1, κy = 0.15), Frame([0., 0., 0.85], MicroSwimmers.I3))
    body = ImplicitExcavateBody(el, groove, 50.0) 
    # jakoba_pars = (a = 3.9, b = 2.2, c = 2.2, a_g = 3.9, b_g = 2.2, c_g = 2.2, p_a = 2, p_b = 2, p_c = 2, z_s = 0.85, θ = 0.0, κ_x = 0.1, κ_y = 0.15)

    b = Part(body, 313, 1117)


    # Construct discretised parts
    excavate = MicroSwimmer([
        Part(body, 313, 1117),
        Part(posterior, 31, 117; location=[-3.8, 0.0, 0.2],orientation=rotation_matrix([0.0, 1.0, 0.0], -π/36)),
        Part(anterior, 31, 117; location=[-3.9, 0., 0],orientation=rotation_matrix([0, 1.0, 0.0], -2π/3) )
    ])


    # swimming
    prob = SwimmingTrajectoryProblem(excavate, eps=0.1, t_final=1.0)
    solve_problem!(prob)
    traj = continue_periodic_trajectory(prob.traj, 100)
    animate(traj, excavate)

    # feeding
    move_boundary!(excavate, [0., 0., 5], MicroSwimmers.I3, 0.0)
    prob = ResistanceProblem(excavate, eps=0.1, wall=true)
    solve_problem!(prob)


    u = FluidVelocity(prob)
    fig = stream(u, range(-15, 15, 50), range(0, 15, 50); plane=:xz)
    move_boundary!(excavate, [0., 5., 0.], rotation_matrix([1., 0., 0.], -π/2), 0.0)
    viz!(fig.content[1], excavate)
    save("FluidVelocity Beat Pattern $mult.png",fig)


    ave_vf = TimeAveragedPlanarVelocityField(prob, range(-15, 15, 50), range(0, 15, 50); plane=:xz)
    fig = stream(ave_vf)
    viz!(fig.content[1], excavate)
    save("TimeAveragedPlanarVelocityField BeatPatten $mult.png",fig)
end


