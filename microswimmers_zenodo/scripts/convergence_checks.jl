using MicroSwimmers
using MicroSwimmersPlots
    
Ns = [54, 216, 864, 3456]
Qs = [864, 3456, 13824, 55296, 221184]

translating_rel_errs = zeros(length(Ns), length(Qs))
for i in eachindex(Ns)
    for j in eachindex(Qs)
        body = CellBody(EllipsoidBody(1.0, 1.0, 1.0), Ns[i], Qs[j])
        prob = ResistanceProblem(body,eps=0.02)
        add_rigid_body_motion!(prob, [1., 0., 0.], zeros(3))
        F, T = total_force_and_torque(prob)
        translating_rel_errs[i,j] = abs(6π - F[1]) / 6π
    end
end
@info "translating sphere relative errors" translating_rel_errs

rotation_rel_errs = zeros(length(Ns), length(Qs))
for i in eachindex(Ns)
    for j in eachindex(Qs)
        body = CellBody(EllipsoidBody(1.0, 1.0, 1.0), Ns[i], Qs[j])
        prob = ResistanceProblem(body,eps=0.02)
        add_rigid_body_motion!(prob, zeros(3), [0., 0., 1.])
        F, T = total_force_and_torque(prob)
        rotation_rel_errs[i,j] = abs(8π - T[3]) / 8π
    end
end

@info "rotating sphere relative errors" rotation_rel_errs


## Flagellum checks

RFT_val_perp = 4π/(log(100) + 0.5)
RFT_val_par = 2π/(log(100) - 0.5) 

N_f = [32, 64, 128, 256]
Q_f = [128, 256, 512, 1024]

perp_rel_errors = zeros(length(N_f), length(Q_f))

RFT_val = RFT_val_perp
for i in eachindex(N_f)
    for j in eachindex(Q_f)
        f = Flagellum(PlanarFlagellum(1., 0., 0., 0., 2π, 2π, 2π, 0.0), N_f[i], Q_f[j])
        prob = ResistanceProblem(f, eps=0.011177755511022043)
        add_rigid_body_motion!(prob, [0. , 0., 1.], zeros(3))
        F, T = total_force_and_torque(prob)
        @info "" F T
        perp_rel_errors[i,j] = abs(F[3] - RFT_val_perp) / RFT_val_perp
    end
end

@info "filament (perpendicular translation)" perp_rel_errors

par_rel_errors = zeros(size(perp_rel_errors))

for i in eachindex(N_f)
    for j in eachindex(Q_f)
        f = Flagellum(PlanarFlagellum(1., 0., 0., 0., 2π, 2π, 2π, 0.0), N_f[i], Q_f[j])
        prob = ResistanceProblem(f, eps=0.011177755511022043)
        add_rigid_body_motion!(prob, [1., 0., 0.], zeros(3))
        F, T = total_force_and_torque(prob)
        par_rel_errors[i,j] = abs(F[1] - RFT_val_par) / RFT_val_par
    end
end

@info "filament (parallel translation)" par_rel_errors

# Optimising epsilon 
perp_errs = []
epss = range(0.0001, 0.025, 500)
for eps in epss
    f = Flagellum(PlanarFlagellum(1., 0., 0., 0., 2π, 2π, 2π, 0.0), 127, 517)
    prob = ResistanceProblem(f, eps=eps)
    add_rigid_body_motion!(prob, [0. , 0., 1.], zeros(3))
    F, T = total_force_and_torque(prob)
    # @info "" F T
    push!(perp_errs, abs(F[3] - RFT_val_perp) / RFT_val_perp)
end

par_errs = []
epss = range(0.0001, 0.025, 500)
for eps in epss
    f = Flagellum(PlanarFlagellum(1., 0., 0., 0., 2π, 2π, 2π, 0.0), 127, 517)
    prob = ResistanceProblem(f, eps=eps)
    add_rigid_body_motion!(prob, [1., 0., 0.], zeros(3))
    F, T = total_force_and_torque(prob)
    push!(par_errs, abs(F[1] - RFT_val_par) / RFT_val_par)
end

function oberbeck_prolate(c, a)
    @assert c > a "c must be the semi-major axis"
    e = sqrt(1 - (a/c)^2)
    L = log((1+e)/(1-e))

    XA = (8/3)*e^3 / (-2e + (1+e^2)*L)
    YA = (16/3)*e^3 / (2e + (3e^2-1)*L)
    XC = (4/3)*e^3*(1-e^2) / (2e - (1-e^2)*L)
    YC = (4/3)*e^3*(2-e^2) / (-2e + (1+e^2)*L)

    return XA, YA, XC, YC
end

function oberbeck_resistance(c, a; μ=1.0)
    XA, YA, XC, YC = oberbeck_prolate(c, a)
    R_axial_trans      = 6π * μ * c * XA
    R_transverse_trans = 6π * μ * c * YA
    R_axial_rot        = 8π * μ * c^3 * XC
    R_transverse_rot   = 8π * μ * c^3 * YC
    return [R_axial_trans, R_transverse_trans, R_axial_rot, R_transverse_rot]
end

function oberbeck_table(; num_force=531, num_quad= 50539,ε=1e-5)
    cases = [(2.0, 1.0), (5.0, 1.0)]
    
    results = map(cases) do (c, a)
        ana = oberbeck_resistance(c, a)
        el  = CellBody(EllipsoidBody(c, a, a), num_force, num_quad)
        R   = grand_resistance_matrix(el, eps=ε)
        num = [
            R[1,1],
            (R[2,2] + R[3,3]) / 2,
            R[4,4],
            (R[5,5] + R[6,6]) / 2
        ]
        errs = @. abs(num - ana) / ana * 100
        (analytical=ana, numerical=num, errors=errs)
    end
     
    results
end

res = oberbeck_table()
@info "ellipsoid matrix relative errors" res[1].errors res[2].errors