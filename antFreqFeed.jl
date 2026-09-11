#Imports
using MicroSwimmers
using MicroSwimmersPlots
using GLMakie
using FastGaussQuadrature
using Statistics
include("excavate_body_design.jl")

function velocity_flux_polar_ellip_z0(u, x, y0, z0, a, b; Nr=20, Nθ=20)
    #Ellipsoid R=1
    R=1

    rs_raw, wrs = gausslegendre(Nr)
    θs_raw, wθs = gausslegendre(Nθ)

    # Affine transforms
    rs = 0.5 * R * (rs_raw .+ 1)  # r ∈ [0, R]
    wrs .= 0.5 * R * wrs          # Jacobian for r

    θs = π * (θs_raw .+ 1)        # θ ∈ [0, 2π]
    wθs .= π * wθs                # Jacobian for θ

    total_flux = 0.0
    for (r, wr) in zip(rs, wrs), (θ, wθ) in zip(θs, wθs)
        z = z0 + r * sin(θ) * b
        y = y0 + r * cos(θ) * a
        vel = u([x, y, z])
        total_flux += vel[1] * r * wr * wθ * a * b # extra r from polar area element
    end

    total_flux
end


function feedingEfficiency(rprob, u)
    fluxes = []
    powers = []
    for t in range(0,1,10)[1:end-1]
        update_boundary!(rprob, t)
        solve_problem!(rprob)
        push!(fluxes, velocity_flux_polar_ellip_z0(u, 0, 0, 10., 2.1, 1.1))
        push!(powers, total_power(rprob))
    end     

    flux = mean(fluxes)
    power = mean(powers)
    @info "" flux
    @info "" power

    e = (flux^2) / power

    e
end

#Posterior flagellum
f = PlanarStandingWaveFlagellum{Float64}(10.0, 6.283185307179586, 0.0, [0.15, 0.0, -0.35, 0.0], [-0.3, 0.4, 0.0, -0.3])

posterior = PlanarVanedFlagellum(f, 0.1, 0.6, .7)

anterior = ThreeDimensionalFlagellum(9., 1.0, 1.25, 0.1, 12.5, 0., 1.0, 1.25, 0.1, 12.5, 0., 0., 0.)
# anterior = ThreeDimensionalFlagellum{Float64}(9.0, 1.0, 0.0, 1.16, 14.0, 0.16, 1.0, 0.8, 0.53, 21.0, -0.16, 0.0, 0.3584073464102069)


##excavate body
# jakoba parameters
el = SuperEllipsoid(3.9, 2.2, 2.2)
groove = Posed(SuperEllipsoid(3.9, 2.2, 2.2; κx = 0.1, κy = 0.15), Frame([0., 0., 0.85], MicroSwimmers.I3))
body = ImplicitExcavateBody(el, groove, 50.0) 
# jakoba_pars = (a = 3.9, b = 2.2, c = 2.2, a_g = 3.9, b_g = 2.2, c_g = 2.2, p_a = 2, p_b = 2, p_c = 2, z_s = 0.85, θ = 0.0, κ_x = 0.1, κ_y = 0.15)

#azimuthal Frequency
aziFrq = collect(0:0.1:2.5)

#elevation Frequency
eleFrq = collect(0:0.1:2.5)

#Number of azi and ele ppoints
nazi = length(aziFrq)
nele = length(eleFrq)

#All Effeciencys
eff = zeros(nele,nazi)

#All fluxes
flu = zeros(nele,nazi)

# #For loop to investigate 
for (col, azi) in enumerate(aziFrq)
    anterior.fᵩ = azi
    for (row, elv) in enumerate(eleFrq)
        anterior.f_θ = elv

        #update problem
        anterior_part = Part(anterior; eps = 0.1, location=[-3.9, 0., 0.25],orientation=rotation_matrix([0, 1.0, 0.0], -2π/3))

        excavate = MicroSwimmer([
            Part(body, 313, 16*313, eps=0.01),
            Part(posterior; eps=0.1, location=[-3.7, 0.0, 0.25],orientation=rotation_matrix([0.0, 1.0, 0.0], -π/36)),
            anterior_part,
            ],
            location=[0., 0., 10.0]
        )

        #Initialise swimming problem 
        rprob = ResistanceProblem(excavate, wall=true)


        #Get fluid velocity
        u = FluidVelocity(rprob)

        #calculate feedingEfficiency
        e = feedingEfficiency(rprob, u)

        # #Effeciency
        println("Effeciency: $e")
        #update Effeciency matrix
        eff[row,col] = e
    end
end

fig = Figure()
ax = Axis(fig[1,1],
    xlabel = "Azimuthal Frequency",
    ylabel = "Elevation Frequency",
)


hm = heatmap!(ax,aziFrq,eleFrq,eff')

Colorbar(fig[1,2],hm,label = "Feeding Effeciency")

save("freqFeedingHEAT.png",fig)



