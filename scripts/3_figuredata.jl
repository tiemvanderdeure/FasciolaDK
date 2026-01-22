### Take julia data and prepare data for figures
# probably merge in a main.jl at some point

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