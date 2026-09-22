using FasciolaDK
using Rasters, DataFrames, CSV, Statistics, LinearAlgebra
using CairoMakie, AlgebraOfGraphics
import GeometryOps as GO, GeoInterface as GI, GADM, GeometryBasics
import Rasters: dims
import FasciolaDK: to_pct_label, Label_subplot, plot_denmark!, vsup_colormatrix, vsup_legend, varname_to_label

global output_date = "2026-09-18" # Dates.today()

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
        3 => "Summer",
        4 => "Autumn",
        1 => "Winter",
        2 => "Spring"
    ]
    colorbarlabels = (
        tavg="Annual average temperature (°C)",
        ppt="Annual total precipitation (mm)",
        ollerenshaw="Annual Ollerenshaw index",
        def="Summer Climate Water Deficit (mm)",
        soil="Summer Soil Moisture (mm)",
    )
    colormaps = (
        tavg = :magma, 
        ppt = :blues, 
        ollerenshaw = :viridis, 
        def = :magma,
        soil = :blues, 
    )
    units = [:tavg, :ppt, :ollerenshaw, :def, :soil] .=> 
        ["Anomaly (°C)", "Anomaly (mm)", "Ollerenshaw index", "Anomaly (mm)", "Anomaly (mm)"]

    # AoG specification
    spec = data(yearly_anomalies) * mapping(
               :year => "",
               units,
               col=AlgebraOfGraphics.dims(1) => renamer(collect(FasciolaDK.envvar_display_names[first.(units)])),
               linestyle=:season => (nonnumeric ∘ renamer(season_names))
           ) * visual(Lines)

    figuregrid = draw(
        spec,
        axis=(limits=(2010, 2023, nothing, nothing), xticks=2010:2:2024, xminorgridvisible=true, titlesize=14),
        figure=(title="Seasonal weather anomalies", titlealign=:center, size=(875, 1200), fontsize=12),
        legend=(labelsize=14, titlesize=14, halign=:left, valign=:top, tellwidth=false)
    )

    fig = figuregrid.figure

    for (i, var) in enumerate(first.(units))
        # plot averages
        ax = Axis(fig[2, i])
        hidedecorations!(ax)
        hidespines!(ax)
        pl = FasciolaDK.plot_denmark!(ax, climatenormals_yearly[var]; colormap=colormaps[var])
        Colorbar(fig[3, i], pl;
            label=colorbarlabels[var],
            labelfont=:bold,
            vertical=false
        )
        # Add labels A-B-C to the subplots
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

save("images/figure1.eps", fig1)

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

save("images/figure2.png", fig2; pt_per_unit=20)

#### Posteriors of effects of weather variables ####
# Figure 3 
posterior_draws = CSV.read("posterior_draws_envvars_$output_date.csv", Tables.rowtable)
posterior_draws = broadcast(posterior_draws) do row
    model = (; row.øko, row.chr, row.slagtid, row.time, row.spacetime)
    return (; model, row.y, row.draw, row.value)
end
posterior_da = DimArray(posterior_draws, (Dim{:model}(DimensionalData.Categorical()), :y, :draw))

focus_models = dims(posterior_da, :model)

posteriors_df = CSV.read("model_runs_posteriors$output_date.csv", DataFrame)
mean_waic_diff = DataFrames.combine(groupby(posteriors_df, :y), [:waic_diff, :dic_diff, :cpo_diff] .=> mean)
envs_sorted = filter(!=("1"), sort(mean_waic_diff, :waic_diff_mean).y)

fig3, figsup_posteriors = map([envs_sorted[1:10], envs_sorted[11:end]]) do envs
    focus_posteriors = posterior_da[y = At(envs)]
    labels = ["Organic + abattoir", "... + spatial", "... + herd ID", "... + time", "... + spatiotemporal"]
    modellabels = [m => l for (m,l) in zip(focus_models, labels)]
    envlabels = [v => varname_to_label(v) for v in envs]

    n_env = length(envs)
    n_model = length(focus_models)
    ylim = quantile(vec(abs.(parent(focus_posteriors))), 0.998) * 1.1

    fig = Figure(size = (max(700, 120 * n_env), 650))

    spec = data(focus_posteriors) * visual(Violin; show_median = true) *
        mapping(
            :y => renamer(envlabels) => "",
            :value => "Posterior effect (log-odds)",
            color = :model => renamer(modellabels) => "Model backbone",
            dodge = :model => renamer(modellabels),
        )

    grid = draw!(fig[1, 1], spec,
        scales(Color = (; palette = from_continuous(:viridis)));
        axis = (
            limits = ((0.5,n_env +0.5), (-ylim, ylim)),
            xticklabelrotation = pi / 4,
            xaxisposition = :top,
            xgridvisible = false,
            xticklabelalign = (:left, :bottom),
            xticksvisible = false
        ),
    )
    legend!(fig[1, 2], grid; title = "Model", labelsize = 18, titlesize = 20)

    ax = grid[1,1].axis
    hlines!(ax, 0; color = :black, linewidth = 1)
    for i in 2:2:n_env
        vspan!(ax, i - 0.5, i + 0.5; color = (:black, 0.06))
    end

    return fig
end

save("images/newfig3.png", fig3)
save("images/figure_S2_posterior_effects.png", figsup_posteriors)


###### Two-variable models ##########
# Figure 4
two_var_df = CSV.read("model_runs_posteriors_twovars_$output_date.csv", DataFrame)

backbone3_single = filter(
    row -> row.øko && row.chr && row.slagtid && !row.time && row.spacetime == "space",
    posteriors_df,
)
heat_vars = envs_sorted[1:7]
n_heat = length(heat_vars)

effect_matrix = fill(NaN, n_heat, n_heat)
effect_matrix025 = fill(NaN, n_heat, n_heat)
effect_matrix975 = fill(NaN, n_heat, n_heat)

map(eachrow(two_var_df)) do row
    i = findfirst(==(row.envvar1), heat_vars)
    j = findfirst(==(row.envvar2), heat_vars)
    effect_matrix[j,i] = row.envvar1_mean
    effect_matrix025[j,i] = row.var"envvar1_quant0.025"
    effect_matrix975[j,i] = row.var"envvar1_quant0.975"
    effect_matrix[i,j] = row.envvar2_mean
    effect_matrix025[i,j] = row.var"envvar2_quant0.025"
    effect_matrix975[i,j] = row.var"envvar2_quant0.975"
end
effect_matrix[diagind(effect_matrix)] .=  backbone3_single.envvar_mean[indexin(heat_vars, backbone3_single.y)]
effect_matrix025[diagind(effect_matrix)] .=  backbone3_single.var"envvar_quant0.025"[indexin(heat_vars, backbone3_single.y)]
effect_matrix975[diagind(effect_matrix)] .=  backbone3_single.var"envvar_quant0.975"[indexin(heat_vars, backbone3_single.y)]

heat_labels = varname_to_label.(heat_vars)

plot_data = effect_matrix[1:end,end:-1:1]
effect_lim = maximum(abs, filter(!isnan, plot_data))
xsize, ysize = size(plot_data)

fig, ax, pl = heatmap(plot_data;
    axis = (
        xticks = (1:xsize, heat_labels[1:xsize]), yticks = (ysize:-1:1, heat_labels[1:ysize]),
        xticklabelrotation = pi / 4, xticklabelalign = (:right, :top),
        ylabel = "Effect of ↓", xlabel = "controlling for →",
        xlabelfont = :bold, ylabelfont = :bold,
        #title = "Two-variable models: effect of each variable controlling for another",
        titlesize = 18,
    ),
    figure = (; size = (900, 800)),
    colormap = Reverse(:RdBu), colorrange = (-effect_lim, effect_lim),
    nan_color = :lightgray,
)

diag_boxes = [Rect2f(k - 0.5, (xsize + 1 - k) - 0.5, 1, 1) for k in 1:xsize]
poly!(ax, diag_boxes; color = :transparent, strokecolor = :black, strokewidth = 1)

labels = cor_to_string.(plot_data)
text!(ax, 
    vec(Tuple.(CartesianIndices(plot_data))); text = vec(labels),
    align = (:center, :center))
Colorbar(fig[1, 2], pl, label = "Posterior mean effect (log-odds)")
fig

save(joinpath("images", "figure4_twovar_effect_heatmap.png"), fig)


############ Posterior estimate of prevalence in space and time ############
# Figure 5
import StatsFuns: logistic
logistic(::Missing) = missing # tiny bit of type piracy to make this work
maybe_quantile(x, args...) = any(ismissing, x) ? missing : quantile(x, args...)

# Spatiotemporal random effect
r_effect_raw = RasterSeries(joinpath("/home/tvd/K/FasciolaDK", "posterior_spatiotemporal.tif"), Dim{:draw}(Int)) |> Rasters.combine
iter_sorted = sort(lookup(r_effect_raw, :draw)) # make sure the order is right here - file system re-orders these so 10 is before 2!
r_effect = set(r_effect_raw[draw=At(iter_sorted)], Rasters.Band => Rasters.format(Dim{:year}(2010:2023))) # fix a dimension

# Environmental variables - but of manual renaming
climate = get_terraclimate((:def,))
climatevars = climate[season=3, year=1:End()-1]
env_vars = RasterStack(
    values(set(climatevars, :year => dims(climatevars, :year) .+ 1)); # fix year dimension
    name=Symbol.(string.(keys(climatevars)) .* "_summer_year1") # fix names
)

# normalization - saved from 1_data_wrangling.jl
normalization = CSV.read("vars_normalized.csv", Tables.columntable; select=collect(keys(env_vars)))
env_vars_norm = maplayers(env_vars, normalization) do v, n
    (v .- n[1]) ./ n[2]
end
year_dim = dims(env_vars_norm, :year)
data_vars = merge(env_vars_norm, 
    RasterStack((; year = DimArray(1:14, (year_dim,)), Intercept = DimArray(fill(1,), ()))))

# Read in fixed effects, def_summer_year1+year+intercept
posterior_env_effect = CSV.read("posterior_effect_def_year_$output_date.csv", Tables.columntable)
posterior_effect_ds = DimStack(posterior_env_effect, (:draw,))
posterior_effects_spatial = maplayers(data_vars, posterior_effect_ds[keys(data_vars)]) do r, e
    @d r .* e
end

effects_combined = mapreduce(identity, (x,y) -> (@d x .+ y), layers(posterior_effects_spatial))

posterior_estimate_logscale = r_effect .+ parent(effects_combined)

posterior_estimate = logistic.(posterior_estimate_logscale)
posterior_estimate_mean = dropdims(mean(posterior_estimate; dims=:draw); dims=:draw)
posterior_estimate_025 = rebuild(
    posterior_estimate_mean; 
    data=parent(maybe_quantile.(eachslice(posterior_estimate; dims=otherdims(posterior_estimate, :draw)), 0.025)))
posterior_estimate_975 = rebuild(
    posterior_estimate_mean; 
    data=parent(maybe_quantile.(eachslice(posterior_estimate; dims=otherdims(posterior_estimate, :draw)), 0.975)))

fig3, figs1, figs2 = map((posterior_estimate_mean, posterior_estimate_025, posterior_estimate_975)) do r
    let fig = Figure(size=(800, 900)),
        nrows = 4,
        colorrange = (0, 0.35),
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

save("images/figure3.eps", fig3)
save("images/figures1_prevalence_lower.png", figs1; pt_per_unit=10)
save("images/figures2_prevalence_upper.png", figs2; pt_per_unit=10)

env_effect_mean = dropdims(mean(posterior_effects_spatial.def_summer_year1; dims=:draw); dims=:draw)
random_effect_mean = dropdims(mean(r_effect; dims=:draw); dims=:draw)

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

########## Supplementary Tables

#### Tables with posteriors of INLA models, grouped by backbone
using SummaryTables, DataFrames, Printf
import WriteDocx as W, StatsBase

to_95_interval(m, l, u) = @sprintf "%.2f (%.2f - %.2f)" m l u
to_95_interval(::Missing, l, u) = ""
add_rank(x::AbstractVector; rev) =
    [@sprintf "%.1f (%d)" xi r for (xi, r) in zip(x, StatsBase.ordinalrank(x; rev))]
to_diff_se(d, se) = @sprintf "%.1f (%.1f)" d se
to_diff_se(::Missing, se) = ""

backbone_labels = ["Organic + abattoir", "... + spatial", "... + herd ID", "... + time", "... + spatiotemporal"]

backbone_groups = groupby(posteriors_df, [:øko, :chr, :slagtid, :time, :spacetime])
backbone_keys = [NamedTuple(k) for k in keys(backbone_groups)]
group_order = [findfirst(==(m), backbone_keys) for m in focus_models]

bordered_cells(v, mask) = [Cell(x; border_bottom=b) for (x, b) in zip(v, mask)]

backbone_tables = map(zip(group_order, backbone_labels)) do (gi, label)
    df = sort(backbone_groups[gi], :mlik, rev=true)
    df.envlabel = varname_to_label.(df.y, sep = " - ")

    fixed_vars = df.time[1] ? ("øko", "year", "envvar") : ("øko", "envvar")
    df = transform(
        df,
        (
            Symbol.(var .* "_" .* ["mean", "quant0.025", "quant0.975"])
            =>
                ByRow(to_95_interval) => Symbol("$(var)_formatted")
            for var in fixed_vars
        )...,
    )


    df.waic_diff_formatted = to_diff_se.(df.waic_diff, df.waic_diff_se)
    df.dic_diff_formatted = to_diff_se.(df.dic_diff, df.dic_diff_se)

    # move the baseline row to the bottom, set off by a rule
    df = vcat(filter(:y => !=("1"), df), filter(:y => ==("1"), df))
    mask = falses(nrow(df))
    mask[end-1] = true

    col_labels = ["Environmental variable", "ΔWAIC (SE)", "ΔDIC (SE)", "Environmental effect (95% CI)", "Organic (95% CI)"]
    col_labels2 = ["Backbone model", "WAIC (SE)", "DIC (SE)", " - ", "Organic (95% CI)"]
    columns =  [:envlabel, :waic_diff_formatted, :dic_diff_formatted, :envvar_formatted, :øko_formatted]
    if df.time[1]
        push!(col_labels, "Year (95% CI)")
        push!(col_labels2, "Year (95% CI)")
        push!(columns, :year_formatted)
    end
    cols = [bordered_cells(df[!, x], mask) for x in columns]
    valuecells = hcat(cols...)

    header = Cell.(col_labels; bold=true, border_bottom=true)
    header2 = Cell.(col_labels2; bold=true, border_bottom=true)
    Table(vcat(header', valuecells[1:end-1,:], header2', valuecells[end:end,:]))
end

slug(s) = strip(replace(lowercase(s), r"[^a-z0-9]+" => "_"), '_')
for (label, tbl) in zip(backbone_labels, backbone_tables)
    open(joinpath("images", "table_posteriors_$(slug(label)).tex"), "w") do io
        show(io, MIME"text/latex"(), tbl)
    end
end

# the "main model"
post_slagteri = CSV.read("slagteri_posterior.csv", DataFrame)
post_slagteri.formatted = to_95_interval.(post_slagteri.mean, post_slagteri.var"0.025quant", post_slagteri.var"0.975quant")
post_slagteri.label = string.(post_slagteri.ID)
headercells = Cell.(["Abattoir ID", "Log odds ratio (95% CI)"]; bold=true, border_bottom=true)
datacells = hcat(Cell.(post_slagteri.label), Cell.(post_slagteri.formatted))
table2 = Table(vcat(headercells', datacells))
open(joinpath("images", "posterior_abattoir_random_effects.tex"), "w") do io
    show(io, MIME"text/latex"(), table2)
end


# Simple figure for graphical abstract
prev = filter(x -> x.scope == "Selected cattle", annual_stats).flukes
yrs = filter(x -> x.scope == "Selected cattle", annual_stats).yr

fig, ax, pl = with_theme(theme_minimal()) do
    lines(yrs, prev, color = "black", axis = (ytickformat=to_pct_label,xticks = LinearTicks(5)), linewidth = 3)
end
save("images/prev_over_time_simple.png", fig; pt_per_unit=10)

######### Supplementary Figures ##############

### Supplementary figure - correlation heatmap
preds = CSV.read("/home/tvd/K/FasciolaDK/predictors.csv", DataFrame)
preds = preds[!, occursin.("year", names(preds))]
preds = preds[!, .!in.(names(preds), Ref(["soil_winter_year2", "def_winter_year2"]))]
# order so variables are grouped first by year, then by variable
preds = preds[!, sortperm(getindex.(names(preds), Ref(1:3)))]

cormat = cor(Matrix(preds))
cormat[1:24 .>= collect(1:24)'] .= NaN
cormat_to_plot = rotl90(cormat[1:end-1, 2:end])
fig, ax, pl = heatmap(cormat_to_plot; 
    axis = (
        xticks = (1:23, names(preds)[end:-1:2]), 
        yticks = (1:23, names(preds)[1:end-1]),
        xticklabelrotation = pi/4,
        xgridvisible = false, ygridvisible = false,
        title = "Correlation between weather variables", titlesize = 22
    ),
    figure = (; size = (1200, 1200)),
    colormap = Reverse(:RdBu), colorrange = (-1,1)
)

text!(ax, 
    vec(Tuple.(CartesianIndices(cormat_to_plot))); text = vec(cor_to_string.(cormat_to_plot)),
    align = (:center, :center))
hidespines!(ax, :r); hidespines!(ax, :t)
Colorbar(fig[1,1], pl, tellheight = false, tellwidth = false,
    halign = 1.0, valign = 0.7, height = Relative(0.5),
    flipaxis = false, label = "Correlation coefficient", labelsize = 18)
fig
save(joinpath("images", "correlation_coefficients.png"), fig)
