using FasciolaDK, DataFrames, Statistics, Rasters, Dates, CSV, StatsBase
import DataFrames: combine # to disambiguate from Rasters.combine

years = 2010:2023

#=
climate = get_terraclimate()

function calculate_risk(rs::RasterStack)
    degree_days = max.(rs.tavg .- 10, 0)
    degree_days .* rs.ppt
end
=#
dyr, dyrtilbes, slagt, fund, beskart, bes_chr, besbrugsart, koord = load_registry()
beskoord = innerjoin(koord, bes_chr, on = :CHRNR)
#### Data wrangling ####

## Combine animals, slagt, and fund - including initial filtering
# to include animals slaughtered at 1-2 years
rename!(slagt, "DATO" => "SLAGTDATO")
fund.value .= true
fundwide = unstack(fund, :SLFUNDKODE, :value)
slfundkoder = Symbol.([371, 374, 375, 377, 379, 381])

dyr_slagt = innerjoin(dyr, slagt, on = :DYR_ID)
dyr_slagt_sel = dyr_slagt[
    (Day(365) .< dyr_slagt.SLAGTDATO .- dyr_slagt.FOEDSELSDATO .< Day(730)) .& 
    (dyr_slagt.FOEDSELSDATO .>= Date(2010)),
    [:DYR_ID, :RACE_ID, :FOEDSELSDATO, :SLAGTDATA_ID, :SLAGTBES_ID, :SLAGTDATO]
]

dyr_slagt_fund = leftjoin(dyr_slagt_sel, fundwide, on = :SLAGTDATA_ID)
for c in slfundkoder
    dyr_slagt_fund[!, c] .= Missings.replace(dyr_slagt_fund[!, c], false)
end
#select!(dyr_slagt_fund, Not(:SLAGTDATA_ID, :FOEDSELSDATO))

#### Combine animal data with information about besætninger
# Including trickier filtering to filter out dairy animals

# No dato til means that brugsart is still active
besbrugsart.DATO_TIL .= Missings.replace(besbrugsart.DATO_TIL, Date(2099))
filter!(r -> r.BRUGSART_ID != 68, besbrugsart) # 68 = slagteri

# keep only the DYR_ID we filtered out before
dyr_bes = dyrtilbes[in.(dyrtilbes.DYR_ID, Ref(Set(dyr_slagt_fund.DYR_ID))), :]
dropmissing!(dyr_bes, :DATOTIL)
#filter!(r -> r.SLAGTBES_ID == r.BES_ID || r.DATOTIL - r.DATOFRA > Day(1), dyr_bes)
sort!(dyr_bes, [:DATOFRA, :DATOTIL])
dyr_bes = dyr_bes[!, [:DYR_ID, :BES_ID, :DATOFRA, :DATOTIL]]

## Combine with brugsart to only select animals that have been meat their whole lives
dyr_bes_brugsart = innerjoin(dyr_bes, besbrugsart, on = :BES_ID) # just join on bes_id even though both have timespans
# now filter out the timespans where there is on overlap
dyr_bes_brugsart.fra = max.(dyr_bes_brugsart.DATOFRA, dyr_bes_brugsart.DATO_FRA) # the actual from
dyr_bes_brugsart.til = min.(dyr_bes_brugsart.DATOTIL, dyr_bes_brugsart.DATO_TIL) # the actual to
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
grps_sel = filter(dyr_bes_brugsart_grp) do g
    allequal(g.øko) && # to make sure øko is boolean
    !any(g.mælk) && # we only want meat cows.
        # OBS: the number of animals filtered out by the above is near 0 in 2011 and then
        # sharply increases. Maybe take another look - possibly include 
    !any(g.BRUGSART_ID .> 20) # all of these are tricky ones like hobbydyr, naturpleje, etc.
end

# again make a dataframe where each animal is 1 row.
# 'history' is the besætninger an animal has been at with from/to, stored in a row
historycols = [:BES_ID, :framonth, :tilmonth]
animals_w_history = DataFrames.combine(grps_sel) do x
    sort!(x, [:framonth, :tilmonth])
    (øko = first(x.øko), history = select(x, historycols))
end

# join on data about slagtedato, 
animals_w_data = innerjoin(animals_w_history, dyr_slagt_fund, on = :DYR_ID)

### Filter data and simplify history
# The rationale here is that we want to know where an animal spent the 2nd summer of its life.

### Now 
slagtbes_count = StatsBase.countmap(animals_w_data.SLAGTBES_ID)
slagtbes_ordered = sort(collect(slagtbes_count), by = last, rev = true)
# the 9 biggest butcheries get their own ID, all other ones are 0 
slagtbes_to_slagteri_id = DataFrame(
    SLAGTBES_ID = first.(slagtbes_ordered), 
    slagteri_id = ifelse.(last.(slagtbes_ordered) .> 20_000, 1:length(slagtbes_ordered), 0)
)
n_slagteri = maximum(slagtbes_to_slagteri_id.slagteri_id) # 9

animals_w_slagteri = innerjoin(animals_w_data, slagtbes_to_slagteri_id, on = :SLAGTBES_ID)

function parse_history(row)
    yr = year(row.FOEDSELSDATO) + 1 # this is the year where we assume transmission _might_ take place
    # kind of crude, but OK for now
    # e.g. what if a cow is born in october and slaughtered october 2 yrs later?
    hist_idx = findfirst(eachrow(row.history)) do r
        r.framonth <= (yr, 3) && r.tilmonth >= (yr, 8)
    end 
    bes_id = isnothing(hist_idx) ? 0 : row.history[hist_idx, :BES_ID]
    return (; yr, bes_id)
end

yr_bes = parse_history.(eachrow(animals_w_slagteri))
animals_w_slagteri.yr = getfield.(yr_bes, :yr)
animals_w_slagteri.BES_ID = getfield.(yr_bes, :bes_id)
animals_filtered = animals_w_slagteri[.!iszero.(animals_w_slagteri.BES_ID), Not([:DYR_ID, :SLAGTDATA_ID, :SLAGTBES_ID, :history])]
# finally filter out any bes that have very few animals - to reduce the complexity of the problem
animals_filtered2 = combine(identity, filter(x -> nrow(x) >= 15, groupby(animals_filtered, :BES_ID)))

## Now that we have all the animals we want, gropu them into cohorts
# Each cohort spent the summer at the same farm and 
cohorts_grps = groupby(animals_filtered2, [:BES_ID, :yr, :slagteri_id])
cohorts = combine(cohorts_grps, Symbol(377) => sum => :positive, :øko => first => :øko, nrow => :count)
bes_included = unique(cohorts.BES_ID)

# get the number of animals in each herd for each year
bessize = mapreduce(vcat, years) do yr
    ids = dyrtilbes[dyrtilbes.DATOFRA .<= Date(yr, 6) .&& .!ismissing.(dyrtilbes.DATOTIL) .&& dyrtilbes.DATOTIL .> Date(yr, 6), :BES_ID]
    cm = countmap(ids)
    DataFrame(BES_ID = collect(keys(cm)), bes_size = collect(values(cm)), yr = yr)
end

cohorts = innerjoin(cohorts, bessize, on = [:BES_ID, :yr])

## Combine with climate data
climateraw = get_terraclimate()
climate = maplayers(climateraw[Ti = Date(2000, 12) .. Date(2023, 11)]) do l
    rebuild(
        l; 
        data = reshape(l, size(l)[1:2]..., 12, :), 
        dims = (dims(l, (X,Y))..., Rasters.format(Dim{:month}(1:12)), Rasters.format(Dim{:year}(2001:2023)))
    )
end

bes_lat_lon = DimVector(tuple.(beskoord.LON, beskoord.LAT), Dim{:BES_ID}(beskoord.BES_ID))
bes_lat_lon_included = bes_lat_lon[BES_ID = At(bes_included)]
besclimate = similar(climate, (dims(bes_lat_lon_included)..., dims(climate, (:month, :year))...))
for bes_id in bes_included 
    (x, y) = bes_lat_lon_included[BES_ID = At(bes_id)]
    dst = @view(besclimate[BES_ID = At(bes_id)])
    src = climate[X = Near(x), Y = Near(y)]
    maplayers(copyto!, dst, src)
end

derived_variables = (
    winter_tavg = (mean, 1:3, :tavg),
    winter_tmin = (mean, 1:3, :tmin),
    winter_prec = (mean, 1:3, :ppt),
    spring_prec = (mean, 4:6, :ppt),
    spring_tavg = (mean, 4:6, :tavg),
    spring_aet =  (mean, 4:6, :aet),
    summer_prec = (mean, 7:9, :ppt),
    summer_tavg = (mean, 7:9, :tavg),
)

besclimate_derived = map(derived_variables) do (f, months, var)
   dropdims(f(besclimate[month = At(months)][var]; dims = :month); dims = :month)
end |> RasterStack

climate_derived = map(derived_variables) do (f, months, var)
    dropdims(f(climate[month = At(months)][var]; dims = :month); dims = :month)
end |> RasterStack

besclimate_array = cat(layers(besclimate_derived)...; dims = Dim{:var}(collect(keys(derived_variables))))

idx_no_climatedata = findall(x -> any(ismissing, x), eachslice(besclimate_array; dims = :BES_ID))
bes_no_climatedata = lookup(dims(besclimate_array, :BES_ID))[idx_no_climatedata]
cohorts_w_climatedata = filter(cohorts) do r
    !(r.BES_ID in bes_no_climatedata)
end

# put everything into a single DF
predictorsdf = map(eachrow(cohorts_w_climatedata)) do r
   merge(besclimate_derived[year = At(r.yr), BES_ID = At(r.BES_ID)], r)
end |> DataFrame 
# normalize
for K in keys(besclimate_derived)
    predictorsdf[!, K] = (predictorsdf[!, K] .- mean(predictorsdf[!, K])) ./ std(predictorsdf[!, K])
end
#predictorsdf.bes_size .= log.(predictorsdf.bes_size)
