# Short script to generate supplementary plot for flukicide usage in Denmark, based on data from VetStat
# https://vetstat.fvst.dk/vetstat/
using DataFrames, CSV, Dates
using StringEncodings

# raw file exported from vetstat
vetstat = CSV.File(
    open("data/export-2026-04-15T10_22_00.329Z.csv", enc"ISO-8859-1");
    normalizenames = true, delim = ';', decimal = ',', groupmark = '.') |> DataFrame
dropmissing!(vetstat, [:CHR_nummer, :Aktivt_stof, :Mængde])
# parse Udleveringsdato as Date
vetstat.date = Date.(vetstat.Udleveringsdato, "dd-mm-yyyy")
vetstat.year = year.(vetstat.date)

filter!(x -> x.Aktivt_stof != "Ivermectin" && x.date < Date(2026, 1, 1), vetstat) # Ivermectin is only used together with Closantel
# Fix negative values??
vetstat.Mængde[vetstat.Mængde .< 0] .*= -1
# Fix way apparent comma problems (for some reason especially for Albendazole?)
vetstat.Mængde[vetstat.Mængde .> 5000 .&& vetstat.Aktivt_stof .== "Albendazol"] ./= 1000


yearly_totals = combine(groupby(vetstat, [:Aktivt_stof, :year]), :Mængde => sum => :Mængde)

cattle_kg = 450

doses = Dict(
    "Triclabendazol" => 12 / 1000 * cattle_kg, # or 30 as pour-on
    "Albendazol" => 0.025 * 0.3 * cattle_kg,
    "Closantel" => 0.2 * 0.1 * cattle_kg,
    "Oxyclozanid" => 34 / 1000 * 0.44 * cattle_kg
)

# Convert grams to number of doses for 500kg cattle
yearly_totals.dose_per_500kg = [doses[comp] for comp in yearly_totals.Aktivt_stof]
yearly_totals.num_doses = yearly_totals.Mængde ./ yearly_totals.dose_per_500kg

# Sort for plotting
sort!(yearly_totals, [:year, :Aktivt_stof])

totals = combine(groupby(yearly_totals, :year), :num_doses => sum => :num_doses, :Mængde => sum => :Mængde)
totals.Aktivt_stof .= "Total"
yearly_totals = vcat(yearly_totals, totals; cols= :union)

# Create plots using Makie and AlgebraOfGraphics
using CairoMakie, AlgebraOfGraphics

# Plot 1: Grams per year
spec1 = data(yearly_totals) * mapping(:year, :Mængde, color=:Aktivt_stof => "Active Compound", group=:Aktivt_stof) * (visual(Lines, linewidth=2.5) + visual(Scatter, markersize=7))
fig1 = draw(spec1; 
    axis=(xlabel="Year", ylabel="Total Grams", title="Flukicide Usage Over Time"))
save("flukicide_grams.png", fig1)

# Plot 2: Number of doses per year (normalized to 500kg cattle)
spec2 = data(yearly_totals) * mapping(:year, :num_doses, color=:Aktivt_stof => "Active Compound", group=:Aktivt_stof) * (visual(Lines, linewidth=2.5) + visual(Scatter, markersize=7))
fig2 = draw(spec2; axis=(xlabel="Year", ylabel="Number of Doses (450kg cattle)", title="Flukicide Usage Over Time"))
save("flukicide_doses.png", fig2)
