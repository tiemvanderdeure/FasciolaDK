using FasciolaDK
using Rasters, DataFrames, CSV, Statistics
using CairoMakie, AlgebraOfGraphics
import GeometryOps as GO, GeoInterface as GI, GADM, GeometryBasics
import Rasters: dims
import FasciolaDK: to_pct_label, Label_subplot, plot_denmark!, vsup_colormatrix, vsup_legend

### Some basic functions and data used in multiple figures
dk_munic = GADM.get("DNK"; depth=2) |> DataFrame
dk_reg = GI.getfeature(GADM.get("DNK"; depth=0)) |> first |> GI.geometry
dk_multipoly = GI.convert(GeometryBasics, dk_reg)
polys_by_area = sort(dk_multipoly.polygons, by=GO.area, rev=true)

###### Climate data over time and space
climate = FasciolaDK.get_terraclimate((:tavg, :ppt, :soil, :def))
ollerenshaw = FasciolaDK.get_ollerenshaw()

climate = merge(climate, (; ollerenshaw))
climatenormals = dropdims(mean(climate; dims=:year); dims=:year)
yearly_anomalies = maplayers(climate, climatenormals) do c, n
    mean.(skipmissing.(eachslice(c .- n; dims=Rasters.commondims(dims(c), (:year, :season)))))
end
yearly_anomalies.ollerenshaw .+= mean(skipmissing(ollerenshaw)) # we actually don't want anomalies here

climatenormals_sel = merge(climatenormals, climatenormals[(:soil, :def)][season=3])
climatenormals_yearly = dropdims(mean(climatenormals_sel; dims=:season); dims=:season)
climatenormals_yearly.ppt .*= 4 # to get seasonal totals

fig1 = let
    season_names = [
        1 => "Winter",
        2 => "Spring",
        3 => "Summer",
        4 => "Autumn"
    ]
    names = (tavg="Temperature", ppt="Precipitation", ollerenshaw="Ollerenshaw index",
        soil="Soil moisture", def="Climate water deficit"
    )
    colorbarlabels = (
        tavg="Annual average temperature (°C)",
        ppt="Annual total precipitation (mm)",
        ollerenshaw="Annual Ollerenshaw index",
        soil="Summer Soil Moisture (mm)",
        def="Summer Climate Water Deficit (mm)"
    )

    # AoG specification
    spec = data(yearly_anomalies) * mapping(
               :year => "",
               [:tavg, :ppt, :ollerenshaw, :soil, :def] .=> ["Anomaly (°C)", "Anomaly (mm)", "Ollerenshaw index", "Anomaly (mm)", "Anomaly (mm)"],
               col=AlgebraOfGraphics.dims(1) => renamer(collect(names)),
               linestyle=:season => renamer(season_names)
           ) * visual(Lines)

    figuregrid = draw(
        spec,
        axis=(limits=(2010, 2023, nothing, nothing), xticks=2010:2:2024, xminorgridvisible=true, titlesize=14),
        figure=(title="Seasonal weather anomalies", titlealign=:center, size=(900, 1200), fontsize=12),
        legend=(labelsize=14, titlesize=14, halign=:left, valign=:top, tellwidth=false)
    )

    fig = figuregrid.figure

    for (i, var, cmap) in zip(1:5, keys(names), [:magma, :blues, :viridis, :blues, :magma])
        ax = Axis(fig[2, i])
        hidedecorations!(ax)
        hidespines!(ax)
        pl = FasciolaDK.plot_denmark!(ax, climatenormals_yearly[var]; colormap=cmap)
        Colorbar(fig[3, i], pl;
            label=colorbarlabels[var],
            labelfont=:bold,
            vertical=false
        )
    end
    # Add labels A-B-C to the subplots
    for i in 1:length(names)
        Label_subplot(fig[1, i, Top()], i, fontsize=14)
    end

    # Move contents around -- kind of annoying
    for i in 1:3
        for j in 1:3
            for c in contents(fig[j, i+3])
                fig[j+3, i] = c
            end
        end
        for c in contents(fig[1, i+3, Top()])
            fig[4, i, Top()] = c
        end
        colsize!(fig.layout, i + 3, Fixed(0))
    end

    # Label Danish regions on the map
    ax = Axis(fig[5, 3])
    hidedecorations!(ax)
    hidespines!(ax)
    plot_denmark!(ax, polys_by_area, color=:lightgray, strokecolor=:black, strokewidth=0.5)

    text!(ax, GO.centroid.(polys_by_area[1:4]); text=["Jutland", "Zealand", "", "Funen"],
        fontsize=16, font=:bold, align=(:center, 0.5))
    Label_subplot(fig[5, 3], 6, fontsize=14)

    tightlimits!(ax)

    fig
end

save("images/figure1.png", fig1; pt_per_unit=10)

#### Figure 2: Prevalence of liver disease in Danish cattle
annual_stats = CSV.read(joinpath("data", "fig_1_annual_stats.csv"), DataFrame)
## Data wrangling - Spatial data
municipality_stats = CSV.read(joinpath("data", "fig_1_municipality_stats.csv"), DataFrame)
municipality_stats.geometry = eval.(Meta.parse.(municipality_stats.geometry)) # silly but it works

# Generate color matrix
fig2 = let years = 2010:2023
    vsup_cmap = vsup_colormatrix(;
        cmap=cgrad(:viridis), n_uncertainty=4,
        max_desat=0.7, pow_desat=1.0, max_light=0.7, pow_light=1
    )
    n_edges = [0, 10, 100, 1000, maximum(municipality_stats.n)]

    munic_colors = FasciolaDK.val_u_to_color(
        municipality_stats.flukes, municipality_stats.n, vsup_cmap;
        colorrange=(0, 0.2),
        u_edges=n_edges
    )

    # AoG specification
    spec = AlgebraOfGraphics.data(annual_stats) *
           mapping(
               :yr => "Year",
               [:flukes, :liver] .=> "Prevalence",
               color=AlgebraOfGraphics.dims(1) =>
                   renamer(["Liver fluke", "Other liver disease"]),
               col=:scope
           ) *
           visual(Lines)

    # Draw the figure
    fig = Figure(size=(600, 800))
    # Line plot for disease over time
    gl_lines = GridLayout(fig[1, 1], alignmode=Outside())

    over_time_grid = draw!(gl_lines, spec,
        axis=(; limits=(first(years), last(years), 0, nothing), ytickformat=to_pct_label)
    )
    legend!(gl_lines[1, 1], over_time_grid,
        position=:bottom, halign=0.15, valign=0.95,# orientation = :horizontal, 
        tellheight=false, tellwidth=false
    )

    # Map for liver flukes in space
    ax_map = Axis(
        fig[2, 1], title="Liver fluke prevalence in selected cattle",
    )
    hidedecorations!(ax_map)
    hidespines!(ax_map)
    poly!(ax_map, municipality_stats.geometry, color=munic_colors, strokewidth=0.2)
    vsup_legend(
        fig[2, 1], vsup_cmap;
        rticks=(collect(range(0, 10; length=5)), [string.(n_edges)[1:4]; ">$(n_edges[4])"]),
        thetaticks=([2.5pi / 4, 2pi / 4, 1.5pi / 4], to_pct_label.(range(0, 0.2; length=3))),
    )

    # labels and layout

    for i in 1:2
        Label_subplot(gl_lines[1, i, Top()], i, fontsize=14)
    end
    Label_subplot(fig[2, 1, Top()], 3, fontsize=14)

    rowsize!(fig.layout, 2, Relative(0.6))
    fig
end;

save("images/figure2.png", fig2; pt_per_unit=10)

##### Posterior estimates of the INLA model for prevalence over space and time
import StatsFuns: logistic
logistic(::Missing) = missing # tiny bit of type piracy to make this work

# Random effect
r_effect_raw = RasterSeries(joinpath("/home/tvd/K/FasciolaDK", "posterior_spatiotemporal.tif"), Dim{:iter}(Int)) |> Rasters.combine
iter_sorted = sort(lookup(r_effect_raw, :iter)) # make sure the order is right here - file system re-orders these so 10 is before 2!
r_effect = set(r_effect_raw[iter=At(iter_sorted)], Rasters.Band => Rasters.format(Dim{:year}(2010:2023))) # fix a dimension

# Environmental variables - but of manual renaming
climate = get_terraclimate((:tavg, :ppt))
climatevars = climate[season=3, year=1:End()-1]
env_vars = RasterStack(
    values(set(climatevars; year=dims(climatevars, :year) .+ 1)); # fix year dimension
    name=Symbol.(string.(keys(climatevars)) .* "_summer_lag1") # fix names
)

# normalization - saved from 1_data_wrangling.jl
normalization = CSV.read("vars_normalized.csv", Tables.columntable; select=collect(keys(env_vars)))
env_vars_norm = maplayers(env_vars, normalization) do v, n
    (v .- n[1]) ./ n[2]
end

posterior_env_effect = CSV.read("posterior_effects_envvars_samples2026-03-11.csv", Tables.rowtable)

env_effects = map(posterior_env_effect) do r
    #mapreduce(+, r[keys(env_vars_norm)], layers(env_vars_norm)) do r, e
    mapreduce(+, r[(:tavg_tm1_summer, :ppt_tm1_summer)], layers(env_vars_norm)) do r, e
        r .* e
    end
end

env_effects_combined = RasterSeries(env_effects, dims(r_effect, :iter)) |> Rasters.combine

intercepts = DimVector(getfield.(posterior_env_effect, :Intercept), dims(r_effect, :iter))

posterior_estimate_logscale = @d r_effect .+ env_effects_combined .+ intercepts strict = false

posterior_estimate = logistic.(posterior_estimate_logscale)
posterior_estimate_mean = dropdims(mean(posterior_estimate; dims=:iter); dims=:iter)
maybe_quantile(x, args...) = any(ismissing, x) ? missing : quantile(x, args...)
posterior_estimate_025 = rebuild(
    posterior_estimate_mean; 
    data=parent(maybe_quantile.(eachslice(posterior_estimate; dims=otherdims(posterior_estimate, :iter)), 0.025)))

posterior_estimate_975 = rebuild(
    posterior_estimate_mean; 
    data=parent(maybe_quantile.(eachslice(posterior_estimate; dims=otherdims(posterior_estimate, :iter)), 0.975)))

fig3, figs1, figs2 = map((posterior_estimate_mean, posterior_estimate_025, posterior_estimate_975)) do r
    let fig = Figure(size=(800, 900)),
        nrows = 4,
        colorrange = (0, 0.4),
        colormap = :plasma

        mean_over_time = dropdims(mean(r; dims=:year); dims=:year)

        for (i, t) in enumerate(dims(r, :year))
            # get row and col index from nrows and i
            row = div(i - 1, nrows) + 1
            col = mod1(i, nrows)
            ax = Axis(fig[row, col], title=string(t))
            hidedecorations!(ax)
            hidespines!(ax)
            pl = FasciolaDK.plot_denmark!(ax, r[year=At(t)]; colormap, colorrange)
        end

        ax = Axis(fig[div(14, nrows)+1, mod1(15, nrows)], title="Average")
        hidedecorations!(ax)
        hidespines!(ax)
        FasciolaDK.plot_denmark!(ax, mean_over_time; colormap, colorrange)

        Colorbar(fig[div(15, nrows)+1, mod1(16, nrows)];
            label="Posterior mean prevalence",
            labelfont=:bold,
            vertical=false,
            colormap, colorrange,
            tellheight=false, tellwidth=false,
            width=Relative(0.8),
            tickformat=to_pct_label
        )
        rowgap!(fig.layout, 5)
        rowgap!(fig.layout, 3, -15) # hack because colorbar overreports height
        colgap!(fig.layout, 5)

        fig
    end
end

save("images/figure3.png", fig3; pt_per_unit=10)
save("images/figures1_prevalence_lower.png", figs1; pt_per_unit=10)
save("images/figures2_prevalence_upper.png", figs2; pt_per_unit=10)

env_effect_mean = dropdims(mean(env_effects_combined; dims=:iter); dims=:iter)
random_effect_mean = dropdims(mean(r_effect; dims=:iter); dims=:iter)

figs3, figs4 = map((env_effect_mean, random_effect_mean), (1, 3)) do r, cr
    let fig = Figure(size=(800, 900)),
        nrows = 4,
        colormap = Reverse(:BrBG_10),
        colorrange = (-cr, cr)

        mean_over_time = dropdims(mean(r; dims=:year); dims=:year)

        for (i, t) in enumerate(dims(r, :year))
            # get row and col index from nrows and i
            row = div(i - 1, nrows) + 1
            col = mod1(i, nrows)
            ax = Axis(fig[row, col], title=string(t))
            hidedecorations!(ax)
            hidespines!(ax)
            pl = FasciolaDK.plot_denmark!(ax, r[year=At(t)]; colormap, colorrange)
        end

        ax = Axis(fig[div(14, nrows)+1, mod1(15, nrows)], title="Average")
        hidedecorations!(ax)
        hidespines!(ax)
        FasciolaDK.plot_denmark!(ax, mean_over_time; colormap, colorrange)

        Colorbar(fig[div(15, nrows)+1, mod1(16, nrows)];
            label="Contribution to log-odds (mean)",
            labelfont=:bold,
            vertical=false,
            colormap, colorrange,
            tellheight=false, tellwidth=false,
            width=Relative(0.8),
        )
        rowgap!(fig.layout, 5)
        rowgap!(fig.layout, 3, -15) # hack because colorbar overreports height
        colgap!(fig.layout, 5)

        fig
    end
end

save("images/figureS3_env_effect.png", figs3; pt_per_unit=10);
save("images/figureS4_random_effect.png", figs4; pt_per_unit=10)

#### Tables with posteriors of INLA models
using SummaryTables, DataFrames, Printf, StatsBase, CSV
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

posteriors_df = CSV.read("model_runs_posteriors2026-03-12.csv", DataFrame)
sort!(posteriors_df, :waic, rev=false)

posteriors_df_formatted = transform(
    posteriors_df,
    (
        Symbol.(var .* "_" .* ["mean", "quant0.025", "quant0.975"])
        =>
            ByRow(to_95_interval) => Symbol("$(var)_formatted")
        for var in ("Intercept", "øko", "besid", "slagteriid", "envvar")
    )...,
)

posteriors_df_formatted.mlik_formatted = add_rank(posteriors_df.mlik; rev=true)
posteriors_df_formatted.waic_formatted = add_rank(posteriors_df.waic; rev=false)
posteriors_df_formatted.dic_formatted = add_rank(posteriors_df.dic; rev=false)
posteriors_df_formatted.formula = pseudo_formula.(eachrow(posteriors_df_formatted))
# drop duplicate formulas, keeping the one with the best mlik
posteriors_df_formatted = unique(posteriors_df_formatted, [:formula, :spacetime])

posteriors_df_formatted.include .= false
posteriors_df_formatted.include[1:10] .= true
posteriors_df_formatted[
    findfirst(x -> x.formula == "organic + herd + abattoir" && x.spacetime == "spacetime", eachrow(posteriors_df_formatted)),
    :include] = true

for sp in unique(posteriors_df_formatted.spacetime)
    posteriors_df_formatted.include[findall(posteriors_df_formatted.spacetime .== sp)[1:3]] .= true
end
posteriors_df_formatted[posteriors_df_formatted.formula.=="null", :include] .= true

labels = Cell.(["Formula", "WAIC (rank)", "m. lik. (rank)", "DIC (rank)", "slope (95% CI)"]; bold=true, border_bottom=true)

grps = groupby(posteriors_df_formatted, :spacetime)
cells = mapreduce(vcat, unique(posteriors_df_formatted.spacetime)) do sp
    grp = grps[(spacetime=sp,)]
    df = grp[grp.include, :]
    datacells = hcat(
        Cell.(df.formula),
        Cell.(df.waic_formatted),
        Cell.(df.mlik_formatted),
        Cell.(df.dic_formatted),
        Cell.(df.envvar_formatted)
    )
    header = repeat(
        [Cell("random effect: " * sp, merge=true, italic=true, border_bottom=true)],
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
headercells = Cell.(["Abattoir", "Log odds ratio (95% CI)"]; bold=true, border_bottom=true)
datacells = hcat(Cell.(post_slagteri.label), Cell.(post_slagteri.formatted))
table2 = Table(vcat(headercells', datacells))

#### Models with 2 variables:
# read the candidate models CSV
cand_df = CSV.read("candidate_models_two_vars_2026-04-29.csv", DataFrame)

# add formatted rank columns
cand_df.mlik_formatted = add_rank(cand_df.mlik; rev=true)
cand_df.waic_formatted = add_rank(cand_df.waic; rev=false)
cand_df.dic_formatted = add_rank(cand_df.dic; rev=false)

cand_df_formatted = transform(
    cand_df,
    (
        Symbol.(var .* "_" .* ["mean", "quant0.025", "quant0.975"])
        =>
            ByRow(to_95_interval) => Symbol("$(var)_formatted")
        for var in ("envvar1", "envvar2")
    )...,
)
cand_df_formatted.slope = cand_df_formatted.envvar1_formatted .* " + " .* cand_df_formatted.envvar2_formatted

sort!(cand_df_formatted, :waic, rev=false)

datacells = hcat(
    Cell.(cand_df_formatted.formula),
    Cell.(cand_df_formatted.waic_formatted),
    Cell.(cand_df_formatted.mlik_formatted),
    Cell.(cand_df_formatted.dic_formatted),
    Cell.(cand_df_formatted.slope)
)

labels

table_two_vars = Table(vcat(labels', datacells))


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
        W.Section([
            W.Paragraph([
                W.Run([W.Text("Table 3. Candidate models with two environmental variables")]),
            ]),
            SummaryTables.to_docx(table_two_vars),
        ]),
    ]),
)

W.save(joinpath("images", "tables.docx"), doc)

# Simple figure for graphical abstract
prev = filter(x -> x.scope == "Selected cattle", annual_stats).flukes
yrs = filter(x -> x.scope == "Selected cattle", annual_stats).yr

fig, ax, pl = with_theme(theme_minimal()) do
    lines(yrs, prev, color = "black", axis = (ytickformat=to_pct_label,xticks = LinearTicks(5)), linewidth = 3)
end
save("images/prev_over_time_simple.png", fig; pt_per_unit=10)