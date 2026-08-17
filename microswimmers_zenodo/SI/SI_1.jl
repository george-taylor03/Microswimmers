using BeatPlanes, MicroSwimmers, MicroSwimmersPlots, JLD2
include("../scripts/theme.jl")
set_theme!(theme_dark())
using GLMakie
GLMakie.activate!()

# SI_1.mp4
@load "results/chlamy_sims_final.jld2"
traj = continue_periodic_trajectory(chlamy_sims[8,16].traj, 100)
animate(
    traj, 
    generate_chlamy(chlamy_sims[8,16]);
    bodycolor=:gray50, 
    color=fc,
    n_fade=1500, 
    step=2, # increase to speed up 
    azimuth=5π/4 
)

