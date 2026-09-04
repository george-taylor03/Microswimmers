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
aziWave = collect(-pi/4:pi/32:pi/4)

#elevation curvatures
eleWave = collect(-pi/4:pi/32:pi/4)

#Number of azi and ele ppoints
nazi = length(aziWave)
nele = length(eleWave)

#All Velocities matrix
vels = zeros(nele,nazi)

#Angular velocites
angVels = similar(vels)

#Torsion
tor = similar(vels)

#Polar Direction
pol = similar(vels)

#azimuthal Direction
aziDir = similar(vels)

#Curvature
curv = similar(vels)

anterior = ThreeDimensionalFlagellum(9., 1.0, 1.25, 0.1, 12.5, 0., 1.0, 1.25, 0.1, 12.5, 0., 0., 0.)
# anterior = ThreeDimensionalFlagellum{Float64}(9.0, 1.0, 0.0, 1.16, 14.0, 0.16, 1.0, 0.8, 0.53, 21.0, -0.16, 0.0, 0.3584073464102069)
design(anterior, limits=(-1., 15., -5., 5., -5., 5.))

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
for (col, azi) in enumerate(aziWave)
    anterior.λᵩ = azi
    for (row, elv) in enumerate(eleWave)
        anterior.λ_θ = elv

        #update parameters
        update_boundary!(excavate,0.0)

        # swimming
        solve_problem!(prob)
        traj = continue_periodic_trajectory(prob.traj, 10)
        # animate(traj, excavate)

        #Fit helix to Trajectory
        helix = fit_helix(traj)

        #Get helix quantites
        vels[row,col] = axis_velocity(helix)
        angVels[row,col] = axis_angular_velocity(helix)
        tor[row,col] = torsion(helix)
        pol[row,col] = axis_polar_angle(helix)
        aziDir[row,col] = axis_azimuthal_angle(helix)
        curv[row,col] = curvature(helix)
    end
end


        

#Create Figure
fig = Figure()
ax = Axis(fig[1,1],
    xlabel = "Azimuthal Wavelength",
    ylabel = "Elevation Wavelength",
)


hm = heatmap!(ax,aziWave,eleWave,vels')

Colorbar(fig[1,2],hm,label = L"Velocity\;  v\;(\mu\mathrm{m}/\text{beat})")

save("WavelengthANDVelocityHEAT.png",fig)

fig = Figure()
ax = Axis(fig[1,1],
    xlabel = "Azimuthal Wavelength",
    ylabel = "Elevation Wavelength",
)


hm = heatmap!(ax,aziWave,eleWave,angVels')

Colorbar(fig[1,2],hm,label = L"Angular \; Velocity \; \omega\;\text{(rad/beat)}")

save("WavelengthANDAngVelocityHEAT.png",fig)

fig = Figure()
ax = Axis(fig[1,1],
    xlabel = "Azimuthal Wavelength",
    ylabel = "Elevation Wavelength",
)


hm = heatmap!(ax,aziWave,eleWave,tor')

Colorbar(fig[1,2],hm,label = L"Torsion \; \tau\;(\mu\mathrm{m}^{-1})")

save("WavelengthANDTorsionHEAT.png",fig)

fig = Figure()
ax = Axis(fig[1,1],
    xlabel = "Azimuthal Wavelength",
    ylabel = "Elevation Wavelength",
)


hm = heatmap!(ax,aziWave,eleWave,pol')

Colorbar(fig[1,2],hm,label = L"Polar \; Angle \; \theta\;\text{(rad)}")

save("WavelengthANDPolarHEAT.png",fig)

fig = Figure()
ax = Axis(fig[1,1],
    xlabel = "Azimuthal Wavelength",
    ylabel = "Elevation Wavelength",
)


hm = heatmap!(ax,aziWave,eleWave,aziDir')

Colorbar(fig[1,2],hm,label = L"Azimuthal \; Direction \; \phi\;\text{(rad)}")

save("WavelengthANDAzidirHEAT.png",fig)

fig = Figure()
ax = Axis(fig[1,1],
    xlabel = "Azimuthal Wavelength",
    ylabel = "Elevation Wavelength",
)


hm = heatmap!(ax,aziWave,eleWave,curv')

Colorbar(fig[1,2],hm,label = L"Curvature\;\kappa\;(\mu\mathrm{m}^{-1})")

save("WavelengthANDcurvatureHEAT.png",fig)

