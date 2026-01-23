using FasciolaDK
using Rasters, DataFrames, CSV, Statistics
using CairoMakie, AlgebraOfGraphics
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

function Label_subplot(gp, n; kw...)
    text = "$('A'+(n-1)))"
    Label(
        gp, text, font = :bold, tellheight = false, tellwidth = false, 
        halign = :left, valign = :top; kw...
    )
end


#### Figure 1: Prevalence of liver disease in Danish cattle
annual_stats = CSV.read(joinpath("data", "fig_1_annual_stats.csv"), DataFrame)
## Data wrangling - Spatial data
municipality_stats = CSV.read(joinpath("data", "fig_1_municipality_stats.csv"), DataFrame)
municipality_stats.geometry = eval.(Meta.parse.(municipality_stats.geometry)) # silly but it works

# Generate color matrix
fig1 = let years = 2010:2023
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

    # labels and layout
        
    for i in 1:2 
        Label_subplot(gl_lines[1,i,Top()], i, fontsize = 14)
    end
    Label_subplot(fig[2,1,Top()], 3, fontsize = 14)

    rowsize!(fig.layout, 2, Relative(0.6))
    fig
end;

save("images/figure1.png", fig1; pt_per_unit = 10)

###### Climate data over time and space
climate = get_terraclimate((:tavg, :ppt))
ollerenshaw = FasciolaDK.get_ollerenshaw()

climate = merge(climate, (; ollerenshaw))
climatenormals = dropdims(mean(climate; dims = :year); dims = :year)
yearly_anomalies = maplayers(climate, climatenormals) do c, n
    mean.(skipmissing.(eachslice(c .- n; dims = Rasters.commondims(dims(c), (:year, :season)))))
end
yearly_anomalies.ollerenshaw .+= mean(skipmissing(ollerenshaw)) # we actually don't want anomalies here

climatenormals_yearly = dropdims(mean(climatenormals; dims = :season); dims = :season)
climatenormals_yearly.ppt .*= 4 # to get seasonal totals

fig2 = let season_names = [            
        1 => "Winter", 
        2 => "Spring", 
        3 => "Summer", 
        4 => "Autumn"
    ],
    names = (tavg = "Temperature", ppt = "Precipitation", ollerenshaw = "Ollerenshaw index"),
    colorbarlabels = (
        tavg = "Annual average temperature (°C)", 
        ppt = "Annual total precipitation (mm)",
        ollerenshaw = "Annual Ollerenshaw index"
    )

    # AoG specification
    spec = data(yearly_anomalies) * mapping(
        :year => "",
        [:tavg, :ppt, :ollerenshaw] .=> ["Anomaly (°C)", "Anomaly (mm)", "Ollerenshaw index"],
        col = AlgebraOfGraphics.dims(1) => renamer(collect(names)),
        linestyle = :season => renamer(season_names)
    ) * visual(Lines)

    figuregrid = draw(
        spec, 
        axis = (limits = (2010, 2023, nothing, nothing), xticks = 2010:2:2024, xminorgridvisible = true, titlesize = 14),
        figure = (title ="Seasonal weather anomalies", titlealign = :center, size = (1100, 800), fontsize = 12)
    )

    fig = figuregrid.figure

    for (i, var, cmap) in zip(1:3, keys(names), [:magma, :blues, :viridis])
        ax = Axis(fig[2,i])
        hidedecorations!(ax); hidespines!(ax)
        pl = FasciolaDK.plot_denmark!(ax, climatenormals_yearly[var]; colormap = cmap)
        Colorbar(fig[3,i], pl; 
            label = colorbarlabels[var], 
            labelfont = :bold,
            vertical = false
        )
    end
    for i in 1:3
        Label_subplot(fig[1,i,Top()], i, fontsize = 14)
        Label_subplot(fig[2,i], i+3, fontsize = 14)
    end

    fig
end;
save("images/figure2.png", fig2; pt_per_unit = 10)


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

save("images/figure3.png", fig3; pt_per_unit = 10)

#### Tables with posteriors of INLA models
using SummaryTables, DataFrames, Printf, StatsBase
import WriteDocx as W

to_95_interval(m, l, u) = @sprintf "%.2f (%.2f - %.2f)" m l u
to_95_interval(::Missing, l, u) = ""
add_rank(x::AbstractVector; rev) = 
    [@sprintf "%.1f (%d)" xi r for (xi, r) in zip(x, StatsBase.ordinalrank(x; rev))]
function pseudo_formula(r)
    vars = String[]
    r.øko && push!(vars, "organic")
    r.besid && push!(vars, "herd")
    r.slagtid && push!(vars, "abattoir")
    r.y != "1" && push!(vars, r.y)
    isempty(vars) ? "null" : join(vars, " + ")
end

posteriors_df = vcat(
    CSV.read("model_runs_posteriors2025-10-31.csv", DataFrame),
    CSV.read("model_runs_posteriors2025-12-02.csv", DataFrame)
)

#posteriors_df = CSV.read("model_runs_posteriors2025-10-31.csv", DataFrame)
sort!(posteriors_df, :waic, rev = false)
sort!(posteriors_df, :mlik, rev = true)

posteriors_df_formatted = transform(
    posteriors_df,
    (
        Symbol.(var .* "_" .* ["mean", "quant0.025", "quant0.975"])
        => ByRow(to_95_interval) => Symbol("$(var)_formatted")
        for var in ("Intercept", "øko", "besid", "slagteriid", "envvar")
    )...,
)

posteriors_df_formatted.mlik_formatted = add_rank(posteriors_df.mlik; rev = true)
posteriors_df_formatted.waic_formatted = add_rank(posteriors_df.waic; rev = false)
posteriors_df_formatted.dic_formatted = add_rank(posteriors_df.dic; rev = false)
posteriors_df_formatted.formula = pseudo_formula.(eachrow(posteriors_df_formatted))

filter(r -> r.y == "1", posteriors_df_formatted)

sort!(posteriors_df_formatted, :mlik, rev = true)

posteriors_df_formatted.include .= false
posteriors_df_formatted.include[1:8] .= true
for sp in unique(posteriors_df_formatted.spacetime)
    posteriors_df_formatted.include[findall(posteriors_df_formatted.spacetime .== sp)[1:3]] .= true
end
posteriors_df_formatted[posteriors_df_formatted.formula .== "null", :include] .= true

labels = Cell.(["Formula", "mlik (rank)", "dic (rank)", "waic (rank)", "slope (95% CI)"]; bold = true, border_bottom = true)

grps = groupby(posteriors_df_formatted, :spacetime)
cells = mapreduce(vcat, unique(posteriors_df_formatted.spacetime)) do sp
    grp = grps[(spacetime = sp,)]
    df = grp[grp.include, :]
    datacells = hcat(
        Cell.(df.formula),
        Cell.(df.mlik_formatted),
        Cell.(df.dic_formatted),
        Cell.(df.waic_formatted),
        Cell.(df.envvar_formatted)
    )
    header = repeat(
        [Cell("random effect: " * sp, merge = true, italic = true, border_bottom = true)], 
    size(datacells, 2)
    )

    return vcat(header', datacells)
end

table1 = Table(vcat(labels', cells))

CSV.write("images/posteriors_table.csv", posteriors_df)

# the "main model"
post_slagteri = CSV.read("slagteri_posterior.csv", DataFrame)
post_slagteri.formatted = to_95_interval.(post_slagteri.mean, post_slagteri.var"0.025quant", post_slagteri.var"0.975quant")
post_slagteri.label = ["Others"; string.('A':('A'+9))]
post_slagteri = post_slagteri[[2:end; 1], :]
headercells = Cell.(["Abattoir", "Log odds ratio (95% CI)"]; bold = true, border_bottom = true)
datacells = hcat(Cell.(post_slagteri.label), Cell.(post_slagteri.formatted))
table2 = Table(vcat(headercells', datacells))

mainpost = posteriors_df_formatted[
    findfirst(r -> r.øko && r.besid && r.slagtid && r.y == "1" && r.spacetime == "spacetime", eachrow(posteriors_df_formatted)), 
    :
]
mainpost.øko



doc = W.Document(
    W.Body([
        W.Section([
            W.Paragraph([
                W.Run([W.Text("Table 1")]),
            ]),
            SummaryTables.to_docx(table1),
        ]),
        W.Section([
            W.Paragraph([
                W.Run([W.Text("Table 2")]),
            ]),
            SummaryTables.to_docx(table2),
        ]),
    ]),
)

W.save(joinpath("images", "tables.docx"), doc)


### MWE for a Makie issue
fig = Figure(size = (800, 800))
for i in 1:15
        row = div(i-1, 4) + 1
    col = mod1(i, 4)
    ax = Axis(fig[row, col], title = string(i))
end
cb = Colorbar(
    fig[4,4], colorrange = (0, 1),
    vertical = false,
    tellheight = false, tellwidth = false,
    label = "Colorbar"
)
fig


### Scatter plot of besætninger and prevalences
positivity_by_bes = combine(
    groupby(filter(r -> 2013 < r.yr < 2017, cohorts), :BES_ID), 
    :positive => sum => :positive, :count => sum => :count, :øko => first => :øko
)
bes_stats = DimStack(
    (prevalence = positivity_by_bes.positive ./ positivity_by_bes.count,
    count = positivity_by_bes.count,
    latlon = bes_lat_lon[At(positivity_by_bes.BES_ID)]),
    Dim{:BES_ID}(positivity_by_bes.BES_ID)
)

fig = Figure()
ax = Axis(fig[1,1])
poly!(ax, dk_munic.geom; color = :transparent, strokewidth = 0.5)
scatter!(ax, vec(bes_stats.latlon); color = vec(bes_stats.prevalence), markersize = 3* sqrt.(vec(bes_stats.count)))
fig

##########
dyr_fund = leftjoin(dyr_slagt, fundwide, on = :SLAGTDATA_ID)
for c in slfundkoder
    dyr_fund[!, c] .= Missings.replace(dyr_fund[!, c], false)
end
dyr_2yr = dyr_fund[(dyr_fund.FOEDSELSDATO .+ Year(2)) .< dyr_fund.SLAGTDATO, :]

#### Plot liver disease over time
counts_over_time = combine(groupby(dyr_2yr, :SLAGTDATO)) do g
    (
        n = nrow(g),
        ikter = count(g[!, Symbol(377)]),
        lever = count(any, eachrow(g[!, slfundkoder]))
    )    
end
sort!(counts_over_time, :SLAGTDATO)

dates = Date(2011, 1, 1):Month(1):Date(2023, 9, 11)

pct_over_time = map(dates) do date
    slagtdays = filter(r -> r.SLAGTDATO > date && r.SLAGTDATO < date + Month(3), counts_over_time)
    tot = sum(slagtdays.n)
    (
        ikter = (sum(slagtdays.ikter) / tot),
        lever = sum(slagtdays.lever) / tot
    )
end

ikter_over_time = getindex.(pct_over_time, :ikter)
lever_over_time = getindex.(pct_over_time, :lever)

fig = Figure()
ax = Axis(
    fig[1, 1], 
    ytickformat = to_pct_label, xtickformat = DateLabel(Date(2011, 2, 15), Month(1))
)

band!(ax, 1:153, zeros(length(lever_over_time)), lever_over_time)
band!(ax, 1:153, zeros(length(lever_over_time)), ikter_over_time)
fig


### Climate data
monthly_climate = maplayers(climate) do l
    map(mean ∘ skipmissing, eachslice(l, dims = (:month,:year)))
end

monthly_avgs = maplayers(x -> mean.(eachslice(x, dims = (:month))), monthly_climate) 
monthly_stds = maplayers(x -> std.(eachslice(x, dims = (:month))), monthly_climate) 

tavg_deviation = vec(monthly_climate.tavg .- monthly_avgs.tavg)
deviations = maplayers(monthly_climate, monthly_avgs, monthly_stds) do c, a, std
    vec((c .- a) ./ std)
end


fig = Figure()
ax = Axis(
    fig[1, 1], 
    xtickformat = DateLabel(Date(2010, 1, 1), Month(1))
)
lines!(ax, deviations.ppt)
fig


### Quick graph of different koder by slagteri over time

# calculate prevalences for each month for each slagteri
#animals_w_slagteri.slagtmonth = length.(range.(Date(2011), animals_filtered2.SLAGTDATO; step = Month(1)))
animals_by_slagteri = groupby(animals_filtered2, [:slagteri_id, :yr])
prevalences = DataFrames.combine(animals_by_slagteri, slfundkoder .=> mean .=> slfundkoder, nrow => :count)
sort!(prevalences, :yr)
prevalences_by_slagteri = groupby(prevalences, :slagteri_id)

# Write the figures
for kode in slfundkoder
    fig = Figure()
    ax = Axis(
        fig[1,1],
        title = koder_slfund[kode],
        limits = (2011 ,2023, 0, nothing),#(1 ,156, 0, nothing), 
        ytickformat = to_pct_label, 
        xticks = 2011:2023#(1:12:156, string.(2011:2023)))
    )
    for K in keys(prevalences_by_slagteri)
        lines!(prevalences_by_slagteri[K].yr, prevalences_by_slagteri[K][!, kode]; 
        color = K.slagteri_id, colorrange = (0,9), colormap = :rainbow, label = string(K.slagteri_id))
    end
    axislegend(ax; nbanks = 2)
    fig
    save(joinpath("images", "0904", string(kode) * "by_slagter_over_time.png"), fig)
end


using Rasters, RasterDataSources, ArchGDAL, GLMakie, GeoMakie
t = Raster(WorldClim{BioClim}, 1, res = "2.5m")
t_f = Raster(WorldClim{Future{BioClim, CMIP6, GFDL_ESM4, SSP370}}, :bio1, date = Date(2090), res = "2.5m")