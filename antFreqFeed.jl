#Imports
using MicroSwimmers
using MicroSwimmersPlots
using GLMakie
using FastGaussQuadrature
using Statistics
include("excavate_body_design.jl")

function velocity_flux_polar2(u, x0, y0, z, R; Nr=20, Nθ=20)
    rs_raw, wrs = gausslegendre(Nr)
    θs_raw, wθs = gausslegendre(Nθ)

    # ϕs_raw, wϕs = gausslegendre(Nϕ)

    # Affine transforms
    rs = 0.5 * R * (rs_raw .+ 1)  # r ∈ [0, R]
    wrs .= 0.5 * R * wrs          # Jacobian for r

    θs = π * (θs_raw .+ 1)        # θ ∈ [0, 2π]
    wθs .= π * wθs                # Jacobian for θ

    total_flux = 0.0
    for (r, wr) in zip(rs, wrs), (θ, wθ) in zip(θs, wθs)
        x = x0 + r * cos(θ)
        y = y0 + r * sin(θ)
        vel = u([x, y, z])
        total_flux += vel[3] * r * wr * wθ  # extra r from polar area element
    end

    total_flux
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

#Radius
r = 3.9

#azimuthal Frequency
aziFrq = collect(0:0.25:5)

#elevation Frequency
eleFrq = collect(0:0.25:5)

#Number of azi and ele ppoints
nazi = length(aziFrq)
nele = length(eleFrq)

#All Effeciencys
eff = zeros(nele,nazi)

# #For loop to investigate 
for (col, azi) in enumerate(aziFrq)
    anterior.fᵩ = azi
    for (i, elv) in enumerate(eleFrq)
        anterior.f_θ = elv
        
        #Start from bottom left (azi and ele -pi/4 to start)
        row = length(eleFrq) - i + 1

        anterior_part = Part(anterior, 31, 117; location=[-3.9, 0., 0.25],orientation=rotation_matrix([0, 1.0, 0.0], -2π/3) )

        excavate = MicroSwimmer([
            Part(body, 313, 3117),
            Part(posterior, 31, 117; location=[-3.7, 0.0, 0.25],orientation=rotation_matrix([0.0, 1.0, 0.0], -π/36)),
            anterior_part
        ])
        # animate(excavate)

        #Initialise swimming problem 
        rprob = ResistanceProblem(excavate, eps=0.1, wall=true)

        # Update Problem
        # update_boundary!(excavate,0.0)
        #Get fluid velocity
        u = FluidVelocity(rprob)
        
        fluxes = []
        powers = []
        for t in range(0,1,10)[1:end-1]
            update_boundary!(rprob, t)
            solve_problem!(rprob)
            push!(fluxes, velocity_flux_polar2(u, 0.0, 0.0, 25, r))
            push!(powers, total_power(rprob))
        end     
        flux = mean(fluxes)
        power = mean(powers)

        e = (flux^2)  / power

        # #Effeciency
        println("Effeciency: $e")
        # println("Power}: $power")

        #update Effeciency matrix
        eff[row,col] = e
        
    end
end

fig = Figure()
ax = Axis(fig[1,1],
    xlabel = "Azimuthal Frequency",
    ylabel = "Elevation Frequency",
)


hm = heatmap!(ax,aziFrq,eleFrq,eff)

Colorbar(fig[1,2],hm,label = "Feeding Effeciency")

save("FrequencyFeedingHEAT.png",fig)

