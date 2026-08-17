using CairoMakie
CairoMakie.activate!()
using ColorSchemes

# colors
bc=(:gray50, 0.3)   # body color
fc= :darkorange
fg= :Oranges  # flagella gradient

const APS_FULLWIDTH_IN = 7.0
const PT_PER_IN = 72  # points per inch is fixed
const APS_FULLWIDTH_PT = APS_FULLWIDTH_IN * PT_PER_IN  # 504 pt
const PX_PER_IN = 96
const APS_FULLWIDTH_PX = round(Int, APS_FULLWIDTH_IN * PX_PER_IN)  # 672 px

function aps_theme()
    Theme(
        fontsize = 8,

        Figure = (
            backgroundcolor = :white,
        ),

        Axis = (
            xlabelsize = 9,
            ylabelsize = 9,
            titlesize = 8,
            xticklabelsize = 8,
            yticklabelsize = 8,
            xticksize = 4,
            yticksize = 4,
            xtickwidth = 0.7,
            ytickwidth = 0.7,
            spinewidth = 0.7,
            xgridvisible = false,
            ygridvisible = false,
        ),

        Legend = (
            labelsize = 8,
            framevisible = false,
        ),

        Lines = (
            linewidth = 1.,
        ),

        Scatter = (
            markersize = 7,
        ),

        Colorbar = (
            ticklabelsize = 8,
            tickwidth = 0.7,
            labelsize = 9,
        ),
    )
end

set_theme!(aps_theme())