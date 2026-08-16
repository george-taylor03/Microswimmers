using MicroSwimmers, MicroSwimmersPlots, GLMakie
using LinearAlgebra, StaticArrays
using Accessors

# ---------------------------------------------------------------------------
# Shapes (unchanged)
# ---------------------------------------------------------------------------

struct SuperEllipsoid{T <: Number} <: ImplicitBodyModel
    semi::SVector{3,T}   # a, b, c
    p::SVector{3,T}      # exponents
    κ::SVector{2,T}      # bend coefficients (x, y) vs z²
end

function SuperEllipsoid(a::T, b::T, c::T; p1::T=2.0, p2::T=2.0, p3::T=2.0, κx::T=0.0, κy::T=0.0) where {T <: Number}
    SuperEllipsoid{T}(SVector{3,T}(a, b, c), SVector{3,T}(p1, p2, p3), SVector{2,T}(κx, κy))
end

param_specs(::Type{<:Vane}) = [
    ParamSpec(:s_start, L"s_0", 0:0.01:1;  group=:vane),
    ParamSpec(:s_end,   L"s_1", 0:0.01:1;  group=:vane),
    ParamSpec(:H,       L"H",   0:0.1:10;  group=:vane),
]

MicroSwimmersPlots.param_specs(::Type{<:SuperEllipsoid}) =  [
        comp(L"a",        0.1:0.05:8,     :body, (@optic _.semi), 1),
        comp(L"b",        0.1:0.05:8,     :body, (@optic _.semi), 2),
        comp(L"c",        0.1:0.05:8,     :body, (@optic _.semi), 3),
        comp(L"\kappa_x", -0.2:0.005:0.2, :body, (@optic _.κ),    1),
        comp(L"\kappa_y", -0.2:0.005:0.2, :body, (@optic _.κ),    2)
]
struct Posed{S <: ImplicitBodyModel} <: ImplicitBodyModel
    shape::S
    frame::Frame
end

Posed(shape; d = zero(SVector{3,Float64}), R = one(SMatrix{3,3,Float64})) =
    Posed(shape, Frame(SVector{3}(d), SMatrix{3,3}(R)))

function implicit_local(s::SuperEllipsoid, y::SVector{3,T}) where {T}
    xc = y[1] + s.κ[1] * y[3]^2
    yc = y[2] + s.κ[2] * y[3]^2
    abs(xc / s.semi[1])^s.p[1] +
    abs(yc / s.semi[2])^s.p[2] +
    abs(y[3] / s.semi[3])^s.p[3] - one(T)
end

MicroSwimmers.implicit(ps::Posed, x::SVector{3}) = implicit_local(ps.shape, inv(ps.frame)(x))
MicroSwimmers.implicit(m::SuperEllipsoid, x::SVector{3}) = implicit_local(m, x)
MicroSwimmers.bounding_radius(ps::Posed) = MicroSwimmers.bounding_radius(ps.shape)
MicroSwimmers.bounding_radius(m::SuperEllipsoid) = 1.5norm(m.semi)
MicroSwimmers.seed(ps::Posed) = ps.frame.location
MicroSwimmers.seed(m::SuperEllipsoid) = zero(SVector{3,Float64})

mutable struct ImplicitExcavateBody{T} <: ImplicitBodyModel
    body::SuperEllipsoid{T}              # in the body frame
    groove::Posed{SuperEllipsoid{T}}     # shape + pose relative to body
    k::T
end

MicroSwimmers.implicit(eb::ImplicitExcavateBody, x::SVector{3}) = smooth_max(
    MicroSwimmers.implicit(eb.body, x),
    -MicroSwimmers.implicit(eb.groove, x),
    eb.k,
)
MicroSwimmers.bounding_radius(eb::ImplicitExcavateBody) = maximum([
    MicroSwimmers.bounding_radius(eb.body),
    MicroSwimmers.bounding_radius(eb.groove) + norm(eb.groove.frame.location),
])
MicroSwimmers.seed(eb::ImplicitExcavateBody) = eb.groove.frame.location

# ---------------------------------------------------------------------------
# Rotations
# ---------------------------------------------------------------------------

Rx(α) = @SMatrix [1.0 0 0; 0 cos(α) -sin(α); 0 sin(α) cos(α)]
Ry(β) = @SMatrix [cos(β) 0 sin(β); 0 1.0 0; -sin(β) 0 cos(β)]
Rz(γ) = @SMatrix [cos(γ) -sin(γ) 0; sin(γ) cos(γ) 0; 0 0 1.0]
Rzyx(a) = Rz(a[3]) * Ry(a[2]) * Rx(a[1])

# inverse of Rzyx, used *once* to seed the angle state from the starting R
function zyx_angles(R)
    β = atan(-R[3,1], hypot(R[3,2], R[3,3]))
    SVector(atan(R[3,2], R[3,3]), β, atan(R[2,1], R[1,1]))
end

# ===========================================================================
# GENERIC MACHINERY (candidate for MicroSwimmersPlots)
# Nothing here knows about excavated bodies — it just wires sliders to a flat
# parameter state that is the single source of truth.
# ===========================================================================

# ParamSpec whose get/set act on a flat *state*, targeting the i-th component
# of the SVector reached by optic `o`.
comp(label, range, group, o, i) = ParamSpec(label, range; group,
    get = st -> o(st)[i],
    set = (st, v) -> set(st, o, MicroSwimmersPlots._setcomp(o(st), i, v)))

# Build one SliderGrid per group, each slider updating a field of `state`.
# Because `state` is the source of truth, every callback is identical — no
# per-group special casing.
function slidergrids!(gridpos, state::Observable, specs; titles = nothing)
    for (col, g) in enumerate(unique(s.group for s in specs))
        gspecs = filter(s -> s.group == g, specs)
        rows = [(label = s.label, range = s.range,
                 startvalue = MicroSwimmersPlots.specget(state[], s)) for s in gspecs]
        r = 1
        if titles !== nothing
            Label(gridpos[1, col], titles(g); font = :bold, tellwidth = false)
            r = 2
        end
        sg = SliderGrid(gridpos[r, col], rows...)
        for (slider, s) in zip(sg.sliders, gspecs)
            on(slider.value) do val
                state[] = MicroSwimmersPlots.specset!(state[], s, val)
            end
        end
    end
end
# ===========================================================================
# EXCAVATE-SPECIFIC (stays here / out of the public package)
# ===========================================================================

# model -> flat state (source of truth). Euler angles are stored explicitly,
# so R is only ever decomposed here, once.
excavate_state(m::ImplicitExcavateBody) = (
    groove = (semi = m.groove.shape.semi, p = m.groove.shape.p, κ = m.groove.shape.κ),
    body   = (semi = m.body.semi,         p = m.body.p,         κ = m.body.κ),
    pose   = (g = m.groove.frame.location, α = zyx_angles(m.groove.frame.orientation)),
    k      = m.k,
)

# flat state -> model. This is the "rotation style" reconstruction applied to
# everything: Rzyx(α) sits alongside the SVector assembly, no special status.
function build_excavate(st)
    T = eltype(st.groove.semi)
    groove = Posed(SuperEllipsoid{T}(st.groove.semi, st.groove.p, st.groove.κ),
                   Frame(st.pose.g, Rzyx(st.pose.α)))
    ImplicitExcavateBody(SuperEllipsoid{T}(st.body.semi, st.body.p, st.body.κ), groove, st.k)
end

# reconstruct just the groove shape (for the left pane)
groove_shape(st) = (T = eltype(st.groove.semi);
    SuperEllipsoid{T}(st.groove.semi, st.groove.p, st.groove.κ))

function excavate_specs()
    gsemi = @optic _.groove.semi
    gκ    = @optic _.groove.κ
    bsemi = @optic _.body.semi
    bκ    = @optic _.body.κ
    gd    = @optic _.pose.g
    gα    = @optic _.pose.α

    [ # groove: intrinsic superellipsoid in its own frame
        comp(L"a",        0.1:0.05:8,     :groove, gsemi, 1),
        comp(L"b",        0.1:0.05:8,     :groove, gsemi, 2),
        comp(L"c",        0.1:0.05:8,     :groove, gsemi, 3),
        comp(L"\kappa_x", -0.2:0.005:0.2, :groove, gκ,    1),
        comp(L"\kappa_y", -0.2:0.005:0.2, :groove, gκ,    2),

        # body
        comp(L"a",        0.1:0.05:8,     :body, bsemi, 1),
        comp(L"b",        0.1:0.05:8,     :body, bsemi, 2),
        comp(L"c",        0.1:0.05:8,     :body, bsemi, 3),
        comp(L"\kappa_x", -0.2:0.005:0.2, :body, bκ,    1),
        comp(L"\kappa_y", -0.2:0.005:0.2, :body, bκ,    2),

        # groove pose relative to body — translation and rotation, uniform
        comp(L"g_x",      -2:0.1:2,       :pose, gd, 1),
        comp(L"g_y",      -2:0.1:2,       :pose, gd, 2),
        comp(L"g_z",      -2:0.1:2,       :pose, gd, 3),
        comp(L"\alpha_x", -π:0.05:π,      :pose, gα, 1),
        comp(L"\alpha_y", -π:0.05:π,      :pose, gα, 2),
        comp(L"\alpha_z", -π:0.05:π,      :pose, gα, 3),
    ]
end

function excavate_body_tool(m::ImplicitExcavateBody;
                            specs = excavate_specs(),
                            samples = (100, 100, 100))   # grid res/axis; lower = snappier drags
    state = Observable(excavate_state(m))

    # left pane depends only on the groove sub-state. ignore_equal_values means
    # body/pose edits produce an equal SuperEllipsoid and are suppressed, so the
    # groove mesh only rebuilds when a groove slider actually moves.
    gshape = map(groove_shape, state; ignore_equal_values = true)

    model = @lift build_excavate($state)

    # coalesce rapid drags: only the latest value gets meshed, off the event task
    gshape_in = Makie.async_latest(gshape, 1)
    model_in  = Makie.async_latest(model, 1)

    fig = Figure()
    ax_groove   = LScene(fig[1, 1])
    ax_excavate = LScene(fig[1, 2])

    groove_mesh = @lift let s = $gshape_in
        mesh_from_function(p -> implicit_local(s, SVector{3}(p));
            origin = Vec3((-1.2 .* s.semi)...),
            widths = Vec3(( 2.4 .* s.semi)...),
            samples)
    end
    body_mesh = @lift let mm = $model_in, a = mm.body.semi
        mesh_from_function(x -> MicroSwimmers.implicit(mm, SVector{3}(x));
            origin = Vec3((-1.05 .* a)...),
            widths = Vec3(( 2.1  .* a)...),
            samples)
    end
    mesh!(ax_groove,   groove_mesh)
    mesh!(ax_excavate, body_mesh)

    slidergrids!(fig[2, :], state, specs)

    display(fig)
    model   # derived; @lift build_excavate keeps it in sync with the sliders
end


## For arranging microswimmers

Rx(α) = @SMatrix [1.0 0 0; 0 cos(α) -sin(α); 0 sin(α) cos(α)]
Ry(β) = @SMatrix [cos(β) 0 sin(β); 0 1.0 0; -sin(β) 0 cos(β)]
Rz(γ) = @SMatrix [cos(γ) -sin(γ) 0; sin(γ) cos(γ) 0; 0 0 1.0]
Rzyx(a) = Rz(a[3]) * Ry(a[2]) * Rx(a[1])

# inverse of Rzyx, used *once* to seed the angle state from the starting R
function zyx_angles(R)
    β = atan(-R[3,1], hypot(R[3,2], R[3,3]))
    SVector(atan(R[3,2], R[3,3]), β, atan(R[2,1], R[1,1]))
end

frame_specs(k, group; grange = -3:0.05:3, αrange = -π:0.05:π) = [
    comp(L"g_x",      grange, group, (@optic _[k].g), 1),
    comp(L"g_y",      grange, group, (@optic _[k].g), 2),
    comp(L"g_z",      grange, group, (@optic _[k].g), 3),
    comp(L"\alpha_x", αrange, group, (@optic _[k].α), 1),
    comp(L"\alpha_y", αrange, group, (@optic _[k].α), 2),
    comp(L"\alpha_z", αrange, group, (@optic _[k].α), 3),
]


function arrange!(sw::MicroSwimmer;
                fps    = 30,
                grange = -3:0.05:3,
                αrange = -π:0.05:π,
                limits = (-5., 5., -5., 5., -5., 5.))

    flag_idx = [i for (i, p) in enumerate(sw.parts) if p.model isa FlagellumModel]
    isempty(flag_idx) && error("no FlagellumModel parts to pose")

    # state: one (g, α) per flagellum; α (Euler angles) is the source of truth
    frame_state(fr::Frame) = (g = fr.location, α = zyx_angles(fr.orientation))
    state = Observable([frame_state(sw.parts[i].frame) for i in flag_idx])

    # a column of six pose sliders per flagellum
    specs = reduce(vcat, frame_specs(k, Symbol(:flagellum_, k); grange, αrange)
                         for k in eachindex(flag_idx))
    title_of = Dict(Symbol(:flagellum_, k) =>
                    "flagellum $(flag_idx[k])  ($(nameof(typeof(sw.parts[flag_idx[k]].model))))"
                    for k in eachindex(flag_idx))

    fig = Figure()
    ax  = Axis3(fig[1, 1]; aspect = :data, limits)
    B   = viz!(ax, sw)

    dt = 1 / fps; t = Ref(0.0)

    # write the current frame-state into the (immutable) parts, then refresh.
    # reuse model + disc so the buffers viz! is watching stay live.
    function repose!()
        st = state[]
        for (k, i) in enumerate(flag_idx)
            p = sw.parts[i]
            sw.parts[i] = Part(p.model, p.disc, Frame(st[k].g, Rzyx(st[k].α)))
        end
        update_boundary!(sw, t[]); update_buffer_observable!(B, sw)
    end
    on(_ -> repose!(), state)     # any pose slider re-poses and redraws

    run = Button(fig[2, 1]; label = "Start/Pause", tellwidth = false)
    slidergrids!(fig[3, 1], state, specs) # titles = g -> title_of[g]
    display(fig)

    isrunning = Observable(false)
    on(_ -> (isrunning[] = !isrunning[]), run.clicks)
    on(run.clicks) do _
        @async while isrunning[]
            isopen(fig.scene) || break
            t[] += dt
            update_boundary!(sw, t[]); update_buffer_observable!(B, sw)
            sleep(dt)
        end
    end
    return fig
end