
Statistics.middle(::Missing, x) = missing

function get_terraclimate(vars = (:aet, :tmax, :tmin, :ppt, :def, :soil), destdir = "data/terraclimate_dk.nc"; force = false)
    if force || !isfile(destdir)
        println("Downloading TerraClimate data")
        urls = NamedTuple(var => "http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_$(var)_1958_CurrentYear_GLOBE.nc" for var in vars)

        lazyrs = RasterStack(urls, lazy = true)
        # fix x-y dimensions
        lazyrs = set(lazyrs, :X => Intervals(Center()), :Y => Intervals(Center()))

        dk = GADM.get("DNK"; depth = 0) |> Rasters.GI.getfeature .|> Rasters.GI.geometry |> first
        terraclimate = crop(lazyrs; to = dk, touches = true)[Ti = Where(>=(Date(2000)))] |> read
        terraclimate_dk = mask(terraclimate; with = dk, boundary = :touches)
        write(destdir, terraclimate_dk, deflatelevel = 2; missingval = Inf, force)
    end
    climate = RasterStack(destdir)
    tavg = middle.(climate.tmax, climate.tmin)
    return process_terraclimate(RasterStack((; tavg, climate...)))
end

function process_terraclimate(climateraw)
    # rebuild to have a dimension with quarter/season
    climate = maplayers(climateraw[Ti = Date(2009, 12) .. Date(2023, 11)]) do l
        rebuild(
            l; 
            data = reshape(l, size(l)[1:2]..., 3, 4, :), 
            dims = map(Rasters.format, (dims(l, (X,Y))..., Dim{:month}(1:3), Dim{:season}(1:4), Dim{:year}(2010:2024)))
        )
    end
    # quarterly climate
    return dropdims(mean(climate; dims = :month); dims = :month)
end
