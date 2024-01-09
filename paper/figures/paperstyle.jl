using Plots
using Printf
using Measures
using ColorSchemes
using LaTeXStrings 

cs1 = ColorScheme([ 
    colorant"#29a2c6",
    colorant"#ffcb18",
    colorant"#ff6d31",
    colorant"#ef597b",
    colorant"#9467bd",
    colorant"#73b66b",
])

default(
    dpi = 600,
    size = (400, 200),                  # tip: x~400 to fit in LaTeX double-column       
    fontfamily = "Computer Modern",
    palette = :tableau_superfishel_stone, #cs1

    grid=false,
    showaxis=true,
    minorticks=:true,
    margin=0mm,
    
    lw=2,

    tickfontsize=8,
    annotationfontsize=8,
    guidefontsize=10,
    legendfontsize=10,
    titlefontsize=10,
    title_location=:left,

    fg_legend=:false,
    label=:none,
)
function no_bg_dark!(plt)
    plot!(
        plt,
        foreground_color_guide=:white,
        foreground_color_axis=:white,
        foreground_color_border=:white,
        tickfontcolor=:white,
        legendfontcolor=:white,
        annotationcolor=:white,
        background_color=nothing
    )
end
