using MicroSwimmers

# # A spherical cell body (semi-axes a = b = c = 1 μm)
body = EllipsoidBody(1.0, 1.0, 1.0)

# discretise the body with N = 213 force points and Q=917 quadrature points
body_disc = Part(body, 213, 917)

# Define a planar flagellum beating pattern (tangent-angle model, Gallagher et al. 2018):
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

# Solve the swimming problem for the rigid-body velocity U, angular velocity Ω,
# and the force distribution
prob = SwimmingProblem(ms)
solve_problem!(prob)

U      = get_U(prob)
Ω      = get_Ω(prob)
forces = get_forces(prob)

# An isolated swimmer has zero net-force and net-torque
F, T = total_force_and_torque(prob)