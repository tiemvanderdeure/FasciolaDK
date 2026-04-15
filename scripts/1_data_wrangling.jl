using FasciolaDK, DataFrames, Statistics, Rasters, Dates, CSV, StatsBase
import DataFrames: combine # to disambiguate from Rasters.combine
import Rasters: dims
years = 2010:2023

# Run this once to download the data
# FasciolaDK.download_terraclimate((:pet, :tmin, :tmax, :ppt, :def, :soil), "data/terraclimate_dk.nc"; silent = false)

dyr, dyrtilbes, slagt, fund, beskart, bes_chr, besbrugsart, koord, udstationeringer = load_registry()
beskoord = innerjoin(koord, bes_chr, on = :CHRNR)
n_records = Dict{String, Int}() # to keep track of when records are filtered out!

dyr.male = .!iseven.(dyr.KONK_ID)

#### Data wrangling ####

## Combine animals, slagt, and fund - including initial filtering
# to include animals slaughtered at 1-2 years
rename!(slagt, "DATO" => "SLAGTDATO")
fund.value .= true
fundwide = unstack(fund, :SLFUNDKODE, :value)
slfundkoder = Symbol.([371, 374, 375, 377, 379, 381])

dyr_slagt = leftjoin(innerjoin(dyr, slagt, on = :DYR_ID), fundwide, on = :SLAGTDATA_ID)
for c in slfundkoder
    dyr_slagt[!, c] .= Missings.replace(dyr_slagt[!, c], false)
end
dyr_slagt.liverdisease = vec(any(Array(dyr_slagt[:, filter(!(==(Symbol(377))), slfundkoder)]); dims = 2))
dyr_slagt.flukes = dyr_slagt[:, Symbol(377)]
n_records["all_animals"] =  nrow(dyr_slagt)

## Initial filtering on just birth and slaughter dates
# We only want animals born in spring and slaughtered in autumn of the next year
dyr_slagt_sel = dyr_slagt[
    # Only select animals that are born in spring of one year, and slaughtered in autumn/fall of the next year
    (
        (year.(dyr_slagt.SLAGTDATO) .== 1 .+ year.(dyr_slagt.FOEDSELSDATO) .&& month.(dyr_slagt.SLAGTDATO) .> 9) .|| 
        (year.(dyr_slagt.SLAGTDATO) .== 2 .+ year.(dyr_slagt.FOEDSELSDATO) .&& month.(dyr_slagt.SLAGTDATO) .<= 3)
    ) .&&
    3 .<= month.(dyr_slagt.FOEDSELSDATO) .<= 6,
    [:DYR_ID, :RACE_ID, :FOEDSELSDATO, :SLAGTDATA_ID, :SLAGTBES_ID, :SLAGTDATO, :liverdisease, :flukes]
]
n_records["age_range"] = nrow(dyr_slagt_sel)

# Filter out animals with race_id 1201, 1202, and 1203 as these are milk cattle
dyr_slagt_nodairyraces = dyr_slagt_sel[dyr_slagt_sel.RACE_ID .> 1203, :]
n_records["not_dairy_race"] = nrow(dyr_slagt_nodairyraces)

#select!(dyr_slagt_fund, Not(:SLAGTDATA_ID, :FOEDSELSDATO))

#### Combine animal data with information about besætninger
# Including trickier filtering to filter out dairy animals

# No dato til means that brugsart is still active
besbrugsart.DATO_TIL .= Missings.replace(besbrugsart.DATOTIL, Date(2099))
filter!(r -> r.BRUGSART_ID != 68, besbrugsart) # 68 = slagteri

# keep only the DYR_ID we filtered out before
dyr_bes = dyrtilbes[in.(dyrtilbes.DYR_ID, Ref(Set(dyr_slagt_nodairyraces.DYR_ID))), :]
dropmissing!(dyr_bes, :DATOTIL)
sort!(dyr_bes, [:DATOFRA, :DATOTIL])
dyr_bes = dyr_bes[!, [:DYR_ID, :BES_ID, :DATOFRA, :DATOTIL]]

#### Select meat cattle
# Combine animals with brugsart based on their besætning
dyr_bes_brugsart = innerjoin(dyr_bes, besbrugsart, on = :BES_ID, makeunique = true) # just join on bes_id even though both have timespans
# Besætninger and their burgsart are time-specific. Filter out where time doesn't overlap
dyr_bes_brugsart.fra = max.(dyr_bes_brugsart.DATOFRA, dyr_bes_brugsart.DATOFRA_1) # the actual from
dyr_bes_brugsart.til = ifelse.(
    ismissing.(dyr_bes_brugsart.DATOTIL_1), 
    dyr_bes_brugsart.DATOTIL, 
    min.(dyr_bes_brugsart.DATOTIL, dyr_bes_brugsart.DATOTIL_1)
 ) # the actual to
dyr_bes_brugsart = dyr_bes_brugsart[
    dyr_bes_brugsart.til .> dyr_bes_brugsart.fra, # sort out rows without any overlap in time
    [:DYR_ID, :BES_ID, :fra, :til, :BRUGSART_ID, :øko, :mælk] # select relevant columns
]
# filter out any besætninger where an animal stayed for a very short time
filter!(x -> x.til - x.fra > Day(10), dyr_bes_brugsart)
dyr_bes_brugsart.framonth = yearmonth.(dyr_bes_brugsart.fra) 
dyr_bes_brugsart.tilmonth = yearmonth.(dyr_bes_brugsart.til)

# Group by animal, and select only certain animals
dyr_bes_brugsart_grp = groupby(dyr_bes_brugsart, :DYR_ID)
# Define which animals to include based on the herds they have belonged to history
# The rationale here is that we want to know where an animal spent the 2nd summer of its life.
animals_sel = DataFrames.combine(dyr_bes_brugsart_grp) do g
        yr = year(findmin(g.fra)[1]) + 1 # the year where we assume transmission _might_ take place
    # has this animal been in this herd the entire transmission season?
    idx = findfirst(eachrow(g)) do r
        r.framonth <= (yr, 3) && r.tilmonth >= (yr, 8)
    end
    include = !isnothing(idx) &&
        !(any(g.mælk)) && # animals that have been in a dairy herd at some point have *much* lower prevalence
        !any(>(20), g.BRUGSART_ID) # all of these are tricky ones like hobbydyr, naturpleje, etc.
    if include
        return (øko = g.øko[idx], BES_ID = g.BES_ID[idx], yr)
    else
        return (øko = false, BES_ID = 0, yr)
    end
end
filter!(x -> !iszero(x.BES_ID), animals_sel)
n_records["registered_as_beef"] = nrow(animals_sel)

# Filter out animals that have been sent to a different herd during the relevant time
dyr_udstationeret = Set(udstationeringer.DYR_ID)
filter!(r -> !in(r.DYR_ID, dyr_udstationeret), animals_sel)
n_records["not_grazing_elsewhere"] = nrow(animals_sel)

# join on data about slagtedato, 
animals_w_data = innerjoin(animals_sel, dyr_slagt_nodairyraces, on = :DYR_ID)

### Now add a slagteri id to each animal corresponding to the slaughterhouse they were slaughtered at
slagtbes_count = StatsBase.countmap(dyr_slagt.SLAGTBES_ID)
slagtbes_ordered = sort(collect(slagtbes_count), by = last, rev = true)
# Any butcheries over 100 000 animals slaughtered get their own ID 
slagtbes_to_slagteri_id = DataFrame(
    SLAGTBES_ID = first.(slagtbes_ordered), 
    slagteri_id = ifelse.(last.(slagtbes_ordered) .> 100_000, 1:length(slagtbes_ordered), 0)
)
n_slagteri = maximum(slagtbes_to_slagteri_id.slagteri_id) # 10

animals_w_slagteri = innerjoin(animals_w_data, slagtbes_to_slagteri_id, on = :SLAGTBES_ID)

# finally filter out any bes that have very few animals - to reduce the complexity of the problem
animals_filtered2 = combine(identity, filter(x -> nrow(x) >= 10, groupby(animals_w_slagteri, :BES_ID)))

## Now that we have all the animals we want, group them into cohorts based on bes and yr
cohorts_grps = groupby(animals_filtered2, [:BES_ID, :yr, :slagteri_id])
cohorts = combine(cohorts_grps, 
    :flukes => sum => :positive, 
    :liverdisease => sum => :positive_other, 
    :øko => first => :øko,
    nrow => :count
)
bes_included = unique(cohorts.BES_ID)

## Read in climate data
climateobs = get_terraclimate((:tavg, :ppt, :def, :soil))
ollerenshaw = FasciolaDK.get_ollerenshaw()

climate_expanded = mapreduce(merge, layers(climateobs)) do x
    RasterStack(x; layersfrom = :season, name = string(Rasters.name(x)) .* "_" .* ["winter", "spring", "summer", "autumn"])
end
climate_nolag = merge(climate_expanded, (; ollerenshaw))

climate_all = FasciolaDK.add_lagged_vars(climate_nolag)

seasons_to_exclude = ("winter_lag1", "autumn")
climate = climate_all[filter(x -> !any(y -> endswith(string(x), y), seasons_to_exclude), keys(climate_all))]

# Combine climate and 
bes_lat_lon = DimVector(tuple.(beskoord.LON, beskoord.LAT), Dim{:BES_ID}(beskoord.BES_ID))
bes_lat_lon_included = bes_lat_lon[BES_ID = At(bes_included)]
besclimate = map(bes_lat_lon_included) do (x,y)
    src = climate[X = Near(x), Y = Near(y)]
end |> RasterSeries |> Rasters.combine

idx_no_climatedata = findall(x -> any(ismissing, x), eachslice(first(layers(besclimate)); dims = :BES_ID))
bes_no_climatedata = lookup(dims(besclimate, :BES_ID))[idx_no_climatedata]
cohorts_w_climatedata = filter(cohorts) do r
    !(r.BES_ID in bes_no_climatedata)
end

# put everything into a single DF
predictorsdf = map(eachrow(cohorts_w_climatedata)) do r
   merge(besclimate[year = At(r.yr), BES_ID = At(r.BES_ID)], r)
end |> DataFrame 

# Normalize all the data columns in predictorsdf
# get means and stds, save those and then Normalize
function normalize!(x)
    m = mean(x)
    s = std(x)
    x .= (x .- m) ./ s
    return (m,s)
end
means_stds = map(keys(besclimate)) do var
    m,s = normalize!(predictorsdf[!, Symbol(var)])
end |> NamedTuple{keys(besclimate)}

CSV.write("vars_normalized.csv", (map(collect, means_stds)))
## In-text numbers
topct(x) = string(round(x * 100, digits=2)) * "%"
in_text_numbers = [
    "Number of animals (total): $(nrow(slagt))",
    "Number of animals (filtered): $(sum(predictorsdf.count))",
    "Pct animals included: $(topct(sum(predictorsdf.count) / nrow(slagt)))",
    "Number of cohorts: $(nrow(predictorsdf))",
    "Number of herds: $(length(unique(predictorsdf.BES_ID)))", 
    "Fluke pct all: $(topct(mean(dyr_slagt.flukes)))",
    "Fluke pct included: $(topct(sum(predictorsdf.positive) / sum(predictorsdf.count)))%",
    "Other liver disease all: $(topct(mean(dyr_slagt.liverdisease)))",
    "Other liver disease included: $(topct(sum(predictorsdf.positive_other) / sum(predictorsdf.count)))%"
]

open(joinpath("in_text_numbers.txt"), "w") do io
    for line in in_text_numbers
        println(io, line)
    end
end

n_records["included"] = sum(predictorsdf.count)
n_records["excluded_birth_slaughter_dates"] = n_records["all_animals"] - n_records["age_range"]
n_records["excluded_dairy_races"] = n_records["age_range"] - n_records["not_dairy_race"]
n_records["excluded_non_beef_herds"] = n_records["not_dairy_race"] - n_records["registered_as_beef"]
n_records["excluded_grazing_elsewhere"] = n_records["registered_as_beef"] - n_records["not_grazing_elsewhere"]
n_records["excluded_bes_under_ten_animals"] = n_records["not_grazing_elsewhere"] - n_records["included"]

open("records_filtering_data.tex", "w") do f
    for (key, value) in n_records
        # Convert snake_case to PascalCase for command name
        write(f, "\\newcommand{\\$(replace(key, "_" => ""))}{$value}\n")
    end
end


## Export predictors
predictors_with_lonlat = innerjoin(predictorsdf, beskoord[!, [:BES_ID, :LON, :LAT]], on = :BES_ID)
CSV.write(joinpath(datadir, "..", "predictors.csv"), predictors_with_lonlat)

## Write files for plotting
import GADM, GeometryBasics, GeoInterface as GI, GeometryOps as GO
## Figure 1 - data over time
all_animals = dyr_slagt[:, [:SLAGTDATO, :liverdisease, :flukes]]
all_animals.yr = year.(all_animals.SLAGTDATO .- Month(3))

all_animals_yearly, selected_animals_yearly = map((all_animals, animals_filtered2)) do d
    combine(groupby(d, :yr)) do g
        n = nrow(g)
        (   
            n = n,
            liver = count(g.liverdisease) / n,
            flukes = count(g.flukes) / n
        )
    end
end
all_animals_yearly.scope .= "All cattle"
selected_animals_yearly.scope .= "Selected cattle"
annual_stats = vcat(all_animals_yearly, selected_animals_yearly)
CSV.write(joinpath("data", "fig_1_annual_stats.csv"), annual_stats)

## Figure 1 - data per municipality
# danish municipalities
dk_munic = DataFrame(GADM.get("DNK"; depth = 2))
dk_munic.geometry = GI.convert.(Ref(GeometryBasics), dk_munic.geom)
dk_munic.munic_id = 1:nrow(dk_munic)
munic_key = CSV.read(joinpath("data", "munics.csv"), DataFrame)
dk_municgrps = leftjoin(dk_munic, munic_key, on = [:munic_id, :NAME_2])
dk_grps = combine(groupby(dk_municgrps, :group)) do g
    # either return just the group and geometry, or use Geometryops to union geometries
    geometry = if nrow(g) == 1
        g.geometry[1]
    else 
       mp = reduce((x,y) -> GO.union(x,y; target = GO.GI.MultiPolygonTrait()), g.geometry)
       GI.convert(GeometryBasics, mp)
    end
    name = first(g.NAME_2)
    return (; geometry, name)
end
dk_grps.group .= 1:nrow(dk_grps)
findfirst(x -> contains(x, "Frederikshavn"), dk_grps.name)

# besætninger and their coordinates
besgeom = beskoord[!, [:BES_ID, :LAT, :LON]]
filter!(r -> r.BES_ID in bes_included, besgeom)
besgeom.geometry = tuple.(besgeom.LON, besgeom.LAT)

# find municipality for each cohort
find_first_munic(co; dk_munic = dk_munic) = findfirst(g -> GO.within(co, g), dk_munic.geometry)
besgeom.group = find_first_munic.(besgeom.geometry; dk_munic = dk_grps)
animals_koords = innerjoin(animals_filtered2, besgeom, on = :BES_ID)
# calculate stats for each municipality
municipality_stats = combine(groupby(animals_koords, :group)) do g
    n = nrow(g)
    nbes = nrow(unique(g[!, [:LAT, :LON]]))
    (   
        n = n,
        nbes = nbes,
        liver = count(g.liverdisease) / n,
        flukes = count(g.flukes) / n
    )
end
municipality_stats = rightjoin(municipality_stats, dk_grps[:, [:group, :geometry]], on = :group)
@assert all(>(5), municipality_stats.nbes)
for c in [:n, :liver, :flukes]
    municipality_stats[!, c] .= Missings.replace(municipality_stats[!, c], 0)
end
CSV.write(joinpath("data", "fig_1_municipality_stats.csv"), municipality_stats)



