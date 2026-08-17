include("excavate_body_design.jl")

## Design posterior flagellum
f = PlanarFlagellum(10., 0., 0., 0.5, π/2, 2π, 2π, 0.0)
design(f, limits=(-1., 15., -5., 5., -5., 5.))

## PlanerstandingwaveFlagellum (Length, Beat Freq, Static curvature, Real part vec amp, imaginary part vec amp, standing wave start / stop, height of wave)
# f = PlanarStandingWaveFlagellum(10., 2π, 0.0, [0., 0.4, 0.25, 0.], [0.0, -0.15, 0.4, 0.0])
f = PlanarStandingWaveFlagellum{Float64}(10.0, 6.283185307179586 * 5, 0.0, [0.15, 0.0, -0.35, 0.0], [-0.3, 0.4, 0.0, -0.3])

posterior = PlanarVanedFlagellum(f, 0.1, 0.6, .7)
design(posterior, limits=(-1., 15., -5., 5., -5., 5.))


## Design 3D anterior flagellum
##threedim (    L::Tfᵩ::T;  Aᵩ::T;  δᵩ::T;  λᵩ::T;  Cᵩ::T          # azimuthal bankf_θ::T; A_θ::T; δ_θ::T; λ_θ::T; C_θ::T          # elevation bankγ::T                                            # overall phaseΔγ::T      )
anterior = ThreeDimensionalFlagellum(9., 1.0, 1.25, 0.1, 12.5, 0., 1.0, 1.25, 0.1, 12.5, 0., 0., 0.)
anterior = ThreeDimensionalFlagellum{Float64}(9.0, 1.0, 0.0, 1.16, 14.0, 0.16, 1.0, 0.8, 0.53, 21.0, -0.16, 0.0, 0.3584073464102069)
design(anterior, limits=(-1., 15., -5., 5., -5., 5.))


## Design excavate body with ventral groove
el = SuperEllipsoid(4., 3., 3.,)
groove = Posed(SuperEllipsoid(4., 3., 3.), Frame([0., 0., 2.0], MicroSwimmers.I3))
body = ImplicitExcavateBody(el, groove, 50.0) # just call it excavate body
excavate_body_tool(body)

# jakoba parameters
el = SuperEllipsoid(3.9, 2.2, 2.2)
groove = Posed(SuperEllipsoid(3.9, 2.2, 2.2; κx = 0.1, κy = 0.15), Frame([0., 0., 0.85], MicroSwimmers.I3))
body = ImplicitExcavateBody(el, groove, 50.0) 
# jakoba_pars = (a = 3.9, b = 2.2, c = 2.2, a_g = 3.9, b_g = 2.2, c_g = 2.2, p_a = 2, p_b = 2, p_c = 2, z_s = 0.85, θ = 0.0, κ_x = 0.1, κ_y = 0.15)
excavate_body_tool(body)

b = Part(body, 313, 1117)
scatter(b.disc.quad_pts, markersize=7)


## Construct discretised parts
excavate = MicroSwimmer([
    Part(body, 313, 3117),
    Part(posterior, 31, 117; location=[-3.7, 0.0, 0.25],orientation=rotation_matrix([0.0, 1.0, 0.0], -π/36)),
    Part(anterior, 31, 117; location=[-3.9, 0., 0.25],orientation=rotation_matrix([0, 1.0, 0.0], -2π/3) )
])
arrange!(excavate, grange=-7:0.05:7)
animate(excavate)
## Arrange flagella

# swimming
prob = SwimmingTrajectoryProblem(excavate, eps=0.1, t_final=1.0, saveat=0.01)
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
save("FluidVelocityotherpara.png",fig)


ave_vf = TimeAveragedPlanarVelocityField(prob, range(-15, 15, 50), range(0, 15, 50); plane=:xz)
fig = stream(ave_vf)
viz!(fig.content[1], excavate)
save("TimeAveragedPlanarVelocityFieldotherpara.png",fig)
