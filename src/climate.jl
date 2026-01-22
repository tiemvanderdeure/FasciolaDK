
_middle(x, y) = Statistics.middle(x,y)
_middle(::Missing, x) = missing

function get_terraclimate(args...; kw...)
    climate = get_terraclimate_raw(args...; kw...)
    return process_terraclimate(climate)
end

function download_terraclimate(vars, destdir = "data/terraclimate_dk.nc"; silent = true)
    silent || println("Downloading TerraClimate data")
    urls = NamedTuple(var => "http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_$(var)_1958_CurrentYear_GLOBE.nc" for var in vars)

    lazyrs = RasterStack(urls, lazy = true)
    # fix x-y dimensions
    lazyrs = set(lazyrs, :X => Intervals(Center()), :Y => Intervals(Center()))

    dk = GADM.get("DNK"; depth = 0) |> Rasters.GI.getfeature .|> Rasters.GI.geometry |> first
    terraclimate = crop(lazyrs; to = dk, touches = true)[Ti = Date(2009, 12) .. Date(2023, 11)] |> read
    terraclimate_dk = mask(terraclimate; with = dk, boundary = :touches)
    write(destdir, terraclimate_dk, deflatelevel = 2; missingval = Inf, force = true)
end

function get_terraclimate_raw(vars = (:pet, :tavg, :ppt, :def, :soil), destdir = "data/terraclimate_dk.nc"; force = false)
    uservars = vars
    if :tavg in vars
        vars = filter(!=(:tavg), vars)
        vars = :tmin in vars ? vars : (vars..., :tmin)
        vars = :tmax in vars ? vars : (vars..., :tmax)
    end
    if force || !isfile(destdir)
        download_terraclimate(vars, destdir; silent = false)
    end
    climate = RasterStack(destdir)
    if :tavg in uservars
        tavg = _middle.(climate.tmax, climate.tmin)
        climate = RasterStack((; tavg, climate...))[uservars]
    end
    return climate
end

function process_terraclimate(climateraw)
    # rebuild to have a dimension with quarter/season
    climate = maplayers(climateraw[Ti = Date(2009, 12) .. Date(2023, 11)]) do l
        rebuild(
            l; 
            data = reshape(l, size(l)[1:2]..., 3, 4, :), 
            dims = map(Rasters.format, (dims(l, (X,Y))..., Dim{:month}(1:3), Dim{:season}(1:4), Dim{:year}(2010:2023)))
        )
    end
    # quarterly climate
    maplayers(climate) do l
        f = Rasters.name(l) === :ppt ? sum : mean
        dropdims(f(l; dims = :month); dims = :month)
    end
end

function get_ollerenshaw()
    climateobsraw = FasciolaDK.get_terraclimate_raw((:ppt, :tavg, :pet))
    prec_days = Raster("data/monthly_prec_days.nc")

    ollerenshaw = broadcast(prec_days, climateobsraw) do n, c
        if ismissing(c.tavg)
            return missing
        elseif c.tavg > 10
            return max(n * ((c.ppt-c.pet)/25.4+5),0)
        else
            return 0.0
        end
    end

    ollerenshaw_yearly = broadcast(Rasters.groupby(ollerenshaw[Ti = 2:End()], Ti => year)) do x
        dropdims(sum(x; dims = Ti); dims = Ti)
    end |> RasterSeries |> Rasters.combine
    ollerenshaw_yearly = set(ollerenshaw_yearly, Ti => Dim{:year})
    return ollerenshaw_yearly
end
