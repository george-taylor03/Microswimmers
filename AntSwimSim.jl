#Imports
using MicroSwimmers
using MicroSwimmersPlots
using GLMakie
using FastGaussQuadrature
using Statistics
include("excavate_body_design.jl")

#Posterior flagellum
f = PlanarStandingWaveFlagellum{Float64}(10.0, 6.283185307179586*5, 0.0, [0.15, 0.0, -0.35, 0.0], [-0.3, 0.4, 0.0, -0.3])

posterior = PlanarVanedFlagellum(f, 0.1, 0.6, .7)

##excavate body
# jakoba parameters
el = SuperEllipsoid(3.9, 2.2, 2.2)
groove = Posed(SuperEllipsoid(3.9, 2.2, 2.2; κx = 0.1, κy = 0.15), Frame([0., 0., 0.85], MicroSwimmers.I3))
body = ImplicitExcavateBody(el, groove, 50.0) 
# jakoba_pars = (a = 3.9, b = 2.2, c = 2.2, a_g = 3.9, b_g = 2.2, c_g = 2.2, p_a = 2, p_b = 2, p_c = 2, z_s = 0.85, θ = 0.0, κ_x = 0.1, κ_y = 0.15)
# excavate_body_tool(body)

#azimuthal curvature
azuCurvs = collect(0:pi/16:pi/4)

#elevation curvatures
eleCurvs = collect(0:pi/16:pi/4)

println(eleCurvs)

#All Velocities matrix
vels = zeros(5,5)

anterior = ThreeDimensionalFlagellum(9., 0., 1.25, 0.1, 12.5, C, 1.0, 1.25, 0.1, 12.5, 0., 0., 0.)
# anterior = ThreeDimensionalFlagellum{Float64}(9.0, 1.0, 0.0, 1.16, 14.0, C - i, 1.0, 0.8, 0.53, 21.0, -0.16, 0.0, 0.3584073464102069)

anterior_part = Part(anterior, 31, 117; location=[-3.9, 0., 0],orientation=rotation_matrix([0, 1.0, 0.0], -2π/3) )

excavate = MicroSwimmer([
    Part(body, 313, 3117),
    Part(posterior, 31, 117; location=[-3.7, 0.0, 0.25],orientation=rotation_matrix([0.0, 1.0, 0.0], -π/36)),
    Part(anterior, 31, 117; location=[-3.9, 0., 0.25],orientation=rotation_matrix([0, 1.0, 0.0], -2π/3) )
])
# animate(excavate)

#Initialise swimming problem 
prob = SwimmingTrajectoryProblem(excavate, eps=0.1, t_final=1.0, saveat=0.01)


# #For loop to investigate 
for azi in range(0, stop = pi/4, step = pi/16)
    anterior.Cᵩ = - azi
    for elv in range(0, stop = pi/4, step = pi/16)
        anterior.C_θ = -elv
        
        #update parameters
        update_boundary!(excavate,0.0)

        # swimming
        solve_problem!(prob)
        traj = continue_periodic_trajectory(prob.traj, 10)
        # animate(traj, excavate)

        #Fit helix to Trajectory
        helix = fit_helix(traj)

        #Append Axis velocity
        vels[Int8(elv/(pi/16))+1,Int8(azi/(pi/16))+1] = axis_velocity(helix)
        #Append Current curvature
    end
end


        

#Create Figure
fig = Figure()
ax = Axis(fig[1,1],
    xlabel = "Curvature Angle (rad)",
    ylabel = "Velocity (μm/beat)"
)

lines!(ax,curvs,vels)

ax = Axis(fig[1,1],
    xlabel = "Azimuthal curvature",
    ylabel = "Elevation curvature",
)


hm = heatmap!(ax,azuCurvs,eleCurvs,vels)

Colorbar(fig[1,2],hm,label = "velocity (μm/beat)")

save("curvatureANDVelocityHEAT.png",fig)

