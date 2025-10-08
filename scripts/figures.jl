using CairoMakie, CSV, AlgebraOfGraphics
import GeometryOps as GO, GeoInterface as GI, GADM, GeometryBasics

koder_slfund = CSV.read(joinpath(datadir, "Koder_slagtefund.csv"), DataFrame)
koder_slfund = NamedTuple(Symbol(r.SLFUNDKODE) => r.SLFUNDTEKST for r in eachrow(koder_slfund))
function to_pct_label(x)
    x1 = x .* 100
    xround = round.(Int, x1)
    x2 = all(isapprox.(x1, xround)) ? xround : x1
    return string.(x2) .* "%"
end

#### Figure 1: Prevalence of liver disease in Danish cattle
annual_stats = CSV.read(joinpath("data", "fig_1_annual_stats.csv"), DataFrame)
## Data wrangling - Spatial data
municipality_stats = CSV.read(joinpath("data", "fig_1_municipality_stats.csv"), DataFrame)
municipality_stats.geometry = eval.(Meta.parse.(municipality_stats.geometry)) # silly but it works

# Generate color matrix
vsup_cmap = vsup_colormatrix(; 
    cmap = cgrad(:viridis), n_uncertainty = 4, 
    max_desat = 0.7, pow_desat = 1.0, max_light = 0.7, pow_light = 1
)
n_edges= [0, 10, 100, 1000, maximum(municipality_stats.n)]

munic_colors = FasciolaDK.val_u_to_color(
    municipality_stats.flukes, municipality_stats.n, vsup_cmap; 
    colorrange = (0,0.2), 
    u_edges = n_edges
);

# AoG specification
spec = AlgebraOfGraphics.data(annual_stats) *
    mapping(
        :yr => "Year",
        [:flukes, :liver] .=> "Prevalence",
        color = AlgebraOfGraphics.dims(1) =>
            renamer(["Liver fluke", "Other liver disease"]),
        col = :scope
    ) *
    visual(Lines)

# Draw the figure
fig = Figure(size = (600, 800))
# Line plot for disease over time
gl_lines = GridLayout(fig[1,1], alignmode = Outside())

over_time_grid = draw!(gl_lines, spec, 
    axis = (; limits = (first(years), last(years), 0, nothing), ytickformat = to_pct_label)
)
legend!(gl_lines[1,1], over_time_grid, 
    position = :bottom, halign = 0.15, valign = 0.95,# orientation = :horizontal, 
    tellheight = false, tellwidth = false
)

# Map for liver flukes in space
ax_map = Axis(
    fig[2,1], title = "Liver fluke prevalence in selected cattle",
)
hidedecorations!(ax_map); hidespines!(ax_map)
poly!(ax_map, municipality_stats.geometry, color = munic_colors, strokewidth = 0.2)
vsup_legend(
    fig[2,1], vsup_cmap;
    rticks = (collect(range(0,10; length = 5)), [string.(n_edges)[1:4]; ">$(n_edges[4])"]),
    thetaticks = ([2.5pi/4, 2pi/4, 1.5pi/4], to_pct_label.(range(0,0.2; length = 3))),
)

rowsize!(fig.layout, 2, Relative(0.6))

fig

save("images/figure1.png", fig)
