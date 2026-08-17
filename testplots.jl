using MicroSwimmers
using MicroSwimmersPlots
using GLMakie

# A spherical cell body (semi-axes a = b = c = 1 μm)
body = EllipsoidBody(1.0, 1.0, 1.0)

# discretise the body with N = 213 force points and Q = 917 quadrature points
body_disc = Part(body, 213, 917)

# Define a planar flagellum beating pattern (tangent-angle model, Gallagher et al. ):
#   θ(s,t) = Cs + (R₀ + R₁ sin(ks/L)) cos(ωt - ϕs/L)
flagellum = PlanarFlagellum(
    10.0,  # L:  length (μm)
    0.0,   # C:  static curvature
    0.6,   # R₀: amplitude envelope
    0.5,   # R₁: spatial modulation of amplitude
    π/2,   # k:  envelope wavenumber
    2π,    # ϕ:  travelling-wave wavenumber
    2π,    # ω:  angular frequency
    0.0,   # δ:  overall phase
)

# discretise the flagellum, attached at the edge of the body on the x-axis
flagellum_disc = Part(
    flagellum, 23, 127,
    location=[1.0, 0.0, 0.0],
    orientation=rotation_matrix([1.0, 0.0, 0.0], 0.0),
)

# Assemble the swimmer from its parts
ms = MicroSwimmer([body_disc, flagellum_disc])

# Visualise the swimmer (hold the left mouse button to rotate the camera, right button to translate)
viz(ms)

# Solve the swimming problem for the rigid-body velocity U, angular velocity Ω,
# and the force distribution
prob = SwimmingTrajectoryProblem(ms)
solve_problem!(prob)

# Visualise the time-averaged velocity field around the swimmer
x_points = range(-10, 15, 50)
y_points = range(-10, 10, 50)
ave_vf = TimeAveragedPlanarVelocityField(prob, x_points, y_points; period=1.0, num_t=30)
stream(ave_vf)