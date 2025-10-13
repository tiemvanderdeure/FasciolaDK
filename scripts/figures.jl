using FasciolaDK
using Rasters, DataFrames, CSV, Statistics
using GLMakie, AlgebraOfGraphics
import GeometryOps as GO, GeoInterface as GI, GADM, GeometryBasics
import Rasters: dims

### Some basic functions and data used in multiple figures
dk_munic = GADM.get("DNK"; depth = 2) |> DataFrame

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

###### Climate data over time and space
climatenormals, climateanomalies = get_terraclimate([:tavg, :aet, :ppt])
yearly_anomalies = maplayers(climateanomalies) do l
    map(mean ∘ skipmissing, eachslice(l, dims = (:season, :year)))
end
anomalies_in_space = dropdims(mean(climatenormals; dims = :season); dims = :season)
anomalies_in_space.ppt .*= 12

fig2 = let season_names = [            
        1 => "Winter", 
        2 => "Spring", 
        3 => "Summer", 
        4 => "Autumn"
    ],
    units = (tavg = "°C", ppt = "mm"),
    names = (tavg = "Temperature", ppt = "Precipitation")

    # AoG specification
    spec = data(yearly_anomalies) * mapping(
        :year => "",
        [:tavg, :ppt] .=> ["Anomaly (°C)", "Anomaly (mm)"],
        col = AlgebraOfGraphics.dims(1) => renamer(collect(names)),
        linestyle = :season => renamer(season_names)
    ) * visual(Lines)

    figuregrid = draw(
        spec, 
        axis = (limits = (2010, 2023, nothing, nothing), xticks = 2010:2:2024, xminorgridvisible = true),
        figure = (title ="Seasonal weather anomalies", titlealign = :center, size = (800, 800), fontsize = 10)
    )

    fig = figuregrid.figure

    colorbarlabels = (
        tavg = "Annual average temperature (°C)", 
        ppt = "Annual total precipitation (mm)"
    )

    for (i, var, cmap) in zip(1:2, [:tavg, :ppt], [:magma, :blues])
        ax = Axis(fig[2,i])
        hidedecorations!(ax); hidespines!(ax)
        pl = FasciolaDK.plot_denmark!(ax, anomalies_in_space[var]; colormap = cmap)
        Colorbar(fig[3,i], pl; 
            label = colorbarlabels[var], 
            labelfont = :bold,
            vertical = false
        )
    end

    fig
end;
save("images/figure2.png", fig2)


##### Posterior estimatse
import StatsFuns: logistic
logistic(::Missing) = missing # tiny bit of type piracy to make this work
posterior_mean = Raster(joinpath("/home/tvd/K", "posterior_space_time.nc"))

rs = RasterSeries(joinpath("/home/tvd/K", "posterior_spatiotemporal.tif"), Dim{:iter}(Int)) |> Rasters.combine
rs = set(rs, Rasters.Band => Ti(2010:2023)) # fix a dimension
vals = logistic.(rs) |> skipmissing |> collect

posterior_mean = dropdims(mean(rs; dims = :iter); dims = :iter)
expected_incidence = logistic.(posterior_mean)

mean_over_time = dropdims(mean(expected_incidence; dims = Ti); dims = Ti)

fig3 = let fig = Figure(size = (800, 900)),
    nrows = 4,
    colorrange = (0, 0.4),
    colormap = :plasma

    for (i, t) in enumerate(dims(expected_incidence, Ti))
        # get row and col index from nrows and i
        row = div(i-1, nrows) + 1
        col = mod1(i, nrows)
        ax = Axis(fig[row, col], title = string(t))
        hidedecorations!(ax); hidespines!(ax)
        pl = FasciolaDK.plot_denmark!(ax, expected_incidence[Ti = At(t)]; colormap, colorrange)
    end

    ax = Axis(fig[div(14, nrows)+1, mod1(15, nrows)], title = "Average")
    hidedecorations!(ax); hidespines!(ax)
    FasciolaDK.plot_denmark!(ax, mean_over_time; colormap, colorrange)

    Colorbar(fig[div(15, nrows)+1, mod1(16, nrows)]; 
        label = "Posterior mean", 
        labelfont = :bold,
        vertical = false,
        colormap, colorrange,
        tellheight = false, tellwidth = false,
        width = Relative(0.8)
    )
    rowgap!(fig.layout, 5)
    rowgap!(fig.layout, 3, -15) # hack because colrobar overreports height
    colgap!(fig.layout, 5)

    fig
end

save("images/figure3.png", fig3)
