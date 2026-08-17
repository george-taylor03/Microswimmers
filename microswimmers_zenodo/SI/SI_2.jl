using BeatPlanes, MicroSwimmers, MicroSwimmersPlots, JLD2
include("../scripts/theme.jl")
set_theme!(theme_dark())
using GLMakie
GLMakie.activate!()

# SI_2.mp4
@load "results/ff_sim.jld2"
animate(
    ff_sim.t, 
    ff_sim.traj, 
    generate_filter_feeder(ff_sim), 
    azimuth=5pi/4, 
    step=5, 
    bodycolor=:gray50, 
    color=fc, 
    linewidth=2 
)