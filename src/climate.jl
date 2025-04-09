
Statistics.middle(::Missing, x) = missing

function get_terraclimate(vars = (:aet, :tmax, :tmin, :ppt), destdir = "data/terraclimate_dk.nc"; force = false)
    climate = if force || !isfile(destdir)
        println("Downloading TerraClimate data")
        urls = NamedTuple(var => "http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_$(var)_1958_CurrentYear_GLOBE.nc" for var in vars)

        lazyrs = RasterStack(urls, lazy = true)
        # fix x-y dimensions
        lazyrs = set(lazyrs, :X => Intervals(Center()), :Y => Intervals(Center()))

        dk = GADM.get("DNK"; depth = 0) |> Rasters.GI.getfeature .|> Rasters.GI.geometry |> first
        terraclimate = crop(lazyrs; to = dk, touches = true)[Ti = Where(>=(Date(2000)))] |> read
        terraclimate_dk = mask(terraclimate; with = dk, boundary = :touches)
        write(destdir, terraclimate_dk, deflatelevel = 2; missingval = Inf, force)
    else
        RasterStack(destdir)
    end
    tavg = middle.(climate.tmax, climate.tmin)
    return RasterStack((; tavg, climate...))
end
