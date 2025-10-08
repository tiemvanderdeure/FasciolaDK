using GADM, DataFrames, Makie, GeoMakie
import Colors: HSL, LCHuv

# Helper functions
function desaturate(col::T, amount)::T where T
    col_lch = LCHuv(col)
    return LCHuv(col_lch.l, col_lch.c * (1-amount), col_lch.h)
end

function lighten(col::T, amount)::T where T
    col_hsl = HSL(col)
    return HSL(col_hsl.h, col_hsl.s, 1 - (1 - col_hsl.l) * (1-amount))
end

# Generate a color matrix
function vsup_colormatrix(;
    cmap, n_uncertainty, 
    max_desat, # 0-1, 1 = most desaturated is grey, 0 = no desaturation
    pow_desat, # higher values to desature more slowly
    max_light, # 0-1 1 = most light is white, 0 = no light
    pow_light  # higher values to lighten more slowly
)

    n_levels = [2^(i-1) for i in 1:n_uncertainty]
    max_levels = n_levels[end]

    rel_col = @. (ceil(Int, (1:max_levels)'/(max_levels/n_levels)) - 0.5) / n_levels
    col_shift = @. 1 - ((1:n_uncertainty) - 1)/(n_uncertainty - 1)

    col = getindex.(Ref(cmap), rel_col)
    col_desat = desaturate.(
        col, 
        max_desat.*col_shift.^pow_desat
    )
    return lighten.(
        col_desat,
        max_light.*col_shift.^pow_light
    )
end

function val_u_to_color(
    v, u, colormatrix;
    colorrange = extrema(v), uncertaintyrange = extrema(u),
    v_edges = range(colorrange...; length = size(colormatrix, 2)+1),
    u_edges = range(uncertaintyrange...; length = size(colormatrix, 1)+1)
    )
    n_uncertainty, n_values = size(colormatrix)
    # Find which bin each value is in
    v_bins = clamp.(searchsortedfirst.(Ref(v_edges), v).-1, 1, n_values)

    u_bins = clamp.(searchsortedfirst.(Ref(u_edges), u).-1, 1, n_uncertainty)

    # Get colors
    getindex.(Ref(colormatrix), u_bins, v_bins)
end

function vsup_legend(gp, vsup_cmap; kw...)
    # Legend axis
    ax_legend = PolarAxis(gp;
        thetalimits = (1.5pi/4, 2.5pi/4),
        rlimits = (0,10),
        rticks = (collect(range(0,10; length = 5)), string.(range(0,1; length = 5))),
        thetaticks = ([1.5pi/4, 2pi/4, 2.5pi/4], ["100%", "50%", "0%"]),
        clip = false,
        height= Relative(0.6), width = Relative(0.5), tellwidth = false,tellheight = false,
        halign = 1.0, valign = 0.9, # found these values to be right through trial and error
        kw...
    )
    # Draw the legend - the meshimage trick is needed because surface plot is bugged
    # See https://github.com/MakieOrg/Makie.jl/issues/5235
    mi = meshimage!(
        ax_legend,
        ax_legend.thetalimits[], 0..10,  rotl90(vsup_cmap), shading = NoShading, 
    )
    mi.plots[1].interpolate[] = false
  #  GLMakie.closeall()
end

#=
# Based on multicolor R package
function vsup_color(
    u, v; 
    cmap, n_uncertainty, n_levels = [2^(i-1) for i in 1:n_uncertainty], 
    max_desat, pow_desat, max_light, pow_light
)
    c = 1-u # certainty
    y = clamp(ceil(Int, c * n_uncertainty), 1, n_uncertainty)
    x = clamp(ceil(Int, v * n_levels[y]), 1, n_levels[y])

    rel_col = ((x - 0.5)/n_levels[y] - 0.5/n_levels[end])*n_levels[end]/(n_levels[end] - 1)
    col_shift = 1 - (y - 1)/(n_uncertainty - 1)

    col = Makie.interpolated_getindex(cmap, rel_col)
    return lighten(
        desaturate(
            col, 
            max_desat*col_shift^pow_desat
        ),
        max_light*col_shift^pow_light
    )
end
=#
## Takes VSUP settings and returns a matrix of colors with
# `n_uncertainty` rows and `2^(n_uncertainty-1)` columns