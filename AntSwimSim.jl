#Imports
using MicroSwimmers
using MicroSwimmersPlots
using GLMakie
using FastGaussQuadrature
using Statistics
include("excavate_body_design.jl")

#Posterior flagellum
f = PlanarStandingWaveFlagellum{Float64}(10.0, 6.283185307179586, 0.0, [0.15, 0.0, -0.35, 0.0], [-0.3, 0.4, 0.0, -0.3])

posterior = PlanarVanedFlagellum(f, 0.1, 0.6, .7)

##excavate body
# jakoba parameters
el = SuperEllipsoid(3.9, 2.2, 2.2)
groove = Posed(SuperEllipsoid(3.9, 2.2, 2.2; κx = 0.1, κy = 0.15), Frame([0., 0., 0.85], MicroSwimmers.I3))
body = ImplicitExcavateBody(el, groove, 50.0) 
# jakoba_pars = (a = 3.9, b = 2.2, c = 2.2, a_g = 3.9, b_g = 2.2, c_g = 2.2, p_a = 2, p_b = 2, p_c = 2, z_s = 0.85, θ = 0.0, κ_x = 0.1, κ_y = 0.15)
# excavate_body_tool(body)

#azimuthal curvature
azuCurvs = collect(-pi/4:pi/32:pi/4)

#elevation curvatures
eleCurvs = collect(-pi/4:pi/32:pi/4)

#All Velocities matrix
vels = zeros(17,17)

#Angular velocites
angVels = zeros(17,17)

#Torsion
tor = zeros(17,17)

#Polar Direction
pol = zeros(17,17)

#azimuthal Direction
aziDir = zeros(17,17)

# anterior = ThreeDimensionalFlagellum(9., 1.0, 1.25, 0.1, 12.5, 0., 1.0, 1.25, 0.1, 12.5, 0., 0., 0.)
anterior = ThreeDimensionalFlagellum{Float64}(9.0, 1.0, 0.0, 1.16, 14.0, 0.16, 1.0, 0.8, 0.53, 21.0, -0.16, 0.0, 0.3584073464102069)

anterior_part = Part(anterior, 31, 117; location=[-3.9, 0., 0.25],orientation=rotation_matrix([0, 1.0, 0.0], -2π/3) )

excavate = MicroSwimmer([
    Part(body, 313, 3117),
    Part(posterior, 31, 117; location=[-3.7, 0.0, 0.25],orientation=rotation_matrix([0.0, 1.0, 0.0], -π/36)),
    anterior_part
])
# animate(excavate)

#Initialise swimming problem 
prob = SwimmingTrajectoryProblem(excavate, eps=0.1, t_final=1.0, saveat=0.01)


# #For loop to investigate 
for azi in range(start = -pi/4, stop = pi/4, step = pi/32)
    anterior.Cᵩ = azi
    #Gets correct index (1 indexing not 0)
    col = round(Int,(azi / (pi/32)+9))
    for elv in range(start = -pi/4, stop = pi/4, step = pi/32)
        anterior.C_θ = elv
        
        row = round(Int,(elv / (pi/32)+9))

        #update parameters
        update_boundary!(excavate,0.0)

        # swimming
        solve_problem!(prob)
        traj = continue_periodic_trajectory(prob.traj, 10)
        # animate(traj, excavate)

        #Fit helix to Trajectory
        helix = fit_helix(traj)

        #Append Axis velocity
        vels[row,col] = axis_velocity(helix)

        angVels[row,col] = axis_angular_velocity(helix)

        tor[row,col] = torsion(helix)

        pol[row,col] = axis_polar_angle(helix)

        aziDir[row,col] = axis_azimuthal_angle(helix)
    end
end


        

#Create Figure
fig = Figure()
ax = Axis(fig[1,1],
    xlabel = "Azimuthal curvature",
    ylabel = "Elevation curvature",
)


hm = heatmap!(ax,azuCurvs,eleCurvs,vels)

Colorbar(fig[1,2],hm,label = "velocity (μm/beat)")

save("curvatureANDVelocityHEAT.png",fig)

fig = Figure()
ax = Axis(fig[1,1],
    xlabel = "Azimuthal curvature",
    ylabel = "Elevation curvature",
)


hm = heatmap!(ax,azuCurvs,eleCurvs,angVels)

Colorbar(fig[1,2],hm,label = "Angular Velocity (μm/beat)")

save("curvatureANDAngVelocityHEAT.png",fig)

fig = Figure()
ax = Axis(fig[1,1],
    xlabel = "Azimuthal curvature",
    ylabel = "Elevation curvature",
)


hm = heatmap!(ax,azuCurvs,eleCurvs,tor)

Colorbar(fig[1,2],hm,label = "Torsion")

save("curvatureANDAngTorsionHEAT.png",fig)

fig = Figure()
ax = Axis(fig[1,1],
    xlabel = "Azimuthal curvature",
    ylabel = "Elevation curvature",
)


hm = heatmap!(ax,azuCurvs,eleCurvs,pol)

Colorbar(fig[1,2],hm,label = "Polar Angle")

save("curvatureANDAngPolarHEAT.png",fig)

fig = Figure()
ax = Axis(fig[1,1],
    xlabel = "Azimuthal curvature",
    ylabel = "Elevation curvature",
)


hm = heatmap!(ax,azuCurvs,eleCurvs,aziDir)

Colorbar(fig[1,2],hm,label = "Azimuthal Direction")

save("curvatureANDAngAzidirHEAT.png",fig)