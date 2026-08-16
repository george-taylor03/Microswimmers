#Imports
using MicroSwimmers
using MicroSwimmersPlots
using GLMakie
using FastGaussQuadrature
using Statistics
include("excavate_body_design.jl")

#Posterior flagellum
f = PlanarStandingWaveFlagellum{Float64}(10.0, 6.283185307179586*5, 0.0, [0.15, 0.0, -0.35, 0.0], [-0.3, 0.4, 0.0, -0.3])

posterior = PlanarVanedFlagellum(f, 0.1, 0.9, 1.0)

##excavate body
# jakoba parameters
el = SuperEllipsoid(3.9, 2.2, 2.2)
groove = Posed(SuperEllipsoid(3.9, 2.2, 2.2; κx = 0.1, κy = 0.15), Frame([0., 0., 0.85], MicroSwimmers.I3))
body = ImplicitExcavateBody(el, groove, 50.0) 
# jakoba_pars = (a = 3.9, b = 2.2, c = 2.2, a_g = 3.9, b_g = 2.2, c_g = 2.2, p_a = 2, p_b = 2, p_c = 2, z_s = 0.85, θ = 0.0, κ_x = 0.1, κ_y = 0.15)
# excavate_body_tool(body)

#Static curvature
C = 0.0

#All curvatures
curvs = []

#All Velocities
vels = []

anterior = ThreeDimensionalFlagellum(9., 0., 1.25, 0.1, 12.5, C, 1.0, 1.25, 0.1, 12.5, 0., 0., 0.)
# anterior = ThreeDimensionalFlagellum{Float64}(9.0, 1.0, 0.0, 1.16, 14.0, C - i, 1.0, 0.8, 0.53, 21.0, -0.16, 0.0, 0.3584073464102069)

anterior_part = Part(anterior, 31, 117; location=[-3.9, 0., 0],orientation=rotation_matrix([0, 1.0, 0.0], -2π/3) )

excavate = MicroSwimmer([
    Part(body, 313, 1117),
    Part(posterior, 31, 117; location=[-3.8, 0.0, 0.2],orientation=rotation_matrix([0.0, 1.0, 0.0], -π/36)),
    anterior_part
])
# animate(excavate)

#Initialise swimming problem 
prob = SwimmingTrajectoryProblem(excavate, eps=0.1, t_final=1.0)

# for i in range(0,stop = pi/4, step = pi/8)
anterior.Cᵩ = pi/4

update_boundary!(excavate,0.0)

# swimming
prob = SwimmingTrajectoryProblem(excavate, eps=0.1, t_final=1.0)
solve_problem!(prob)
traj = continue_periodic_trajectory(prob.traj, 10)
# animate(traj, excavate)

#Fit helix to Trajectory
helix = fit_helix(traj)

#Append Axis velocity
push!(vels,axis_velocity(helix))
#Append Current curvature
push!(curvs, C)

#Update anterior

# end

#Great Figure
fig = Figure()
ax = Axis(
    xlabel = "curvature angle",
    ylabel! = "velocity"
)

lines!(ax,curvs,vels)

