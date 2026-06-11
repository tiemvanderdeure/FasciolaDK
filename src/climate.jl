
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
    terraclimate = crop(lazyrs; to = dk, touches = true)[Ti = Date(2008, 12) .. Date(2023, 11)] |> read
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
    climate = maplayers(climateraw[Ti = Date(2008, 12) .. Date(2023, 11)]) do l
        rebuild(
            l; 
            data = reshape(l, size(l)[1:2]..., 3, 4, :), 
            dims = map(Rasters.format, (dims(l, (X,Y))..., Dim{:month}(1:3), Dim{:season}(1:4), Dim{:year}(2009:2023)))
        )
    end
    # quarterly climate
    climate_quarterly = maplayers(climate) do l
        f = Rasters.name(l) === :ppt ? sum : mean
        dropdims(f(l; dims = :month); dims = :month)
    end
#=
    climate_tm1_quarterly = set(climate_quarterly[year = 1:End()-1], year = Dim{:year}(2010:2023))
    climate_tm1_quarterly_rn = RasterStack(NamedTuple(Symbol.(string.(keys(climate_quarterly)) .* "_tm1") .=> values(climate_tm1_quarterly)))

    climate_quarterly_tm0 = climate_quarterly[year = 2:End()]

    return merge(climate_quarterly_tm0, climate_tm1_quarterly_rn)
    =#
end

function add_lagged_vars(climate)
    climate_tm1_quarterly = set(climate[year = 1:End()-1], year = Dim{:year}(val(lookup(climate, :year))[2:end]))
    climate_tm1_quarterly_yr1 = RasterStack(NamedTuple(
        Symbol(string(k) * "_year1") => climate_tm1_quarterly[k] for k in keys(climate)))
    climate_quarterly_yr2 = RasterStack(NamedTuple(Symbol(string(k) * "_year2") => climate[k] for k in keys(climate)))[year = 2:End()]
    return merge(climate_tm1_quarterly_yr1, climate_quarterly_yr2)
end


function get_ollerenshaw(force = false)
    climateobsraw = get_terraclimate_raw((:ppt, :tavg, :pet))[Ti = Where(>=(Date(2008, 12)))]
    prec_days = _get_prec_days(force)

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


##### Rain days for ollerenshaw

import HTTP
function _maybe_download_eobs(force = false)
    filename = "rr_ens_spread_0.1deg_reg_v32.0e.nc"
    uri = HTTP.URI(; scheme = "https", host = "knmi-ecad-assets-prd.s3.amazonaws.com", path = "/ensembles/data/Grid_0.1deg_reg_ensemble/$filename")
    dir = joinpath(ENV["RASTERDATASOURCES_PATH"], "E-OBS")
    isdir(dir) || mkdir(dir)
    filepath = joinpath(dir, filename)
    if !isfile(filepath)
        mkpath(dirname(filepath))
        @info "Starting download for $uri"
        try
            HTTP.download(string(uri), filepath)
        catch e
            # Remove anything that was downloaded before the error
            isfile(filepath) && rm(filepath)
            throw(e)
        end
    end
    return filepath
end


using Rasters.Lookups, Dates
# little trick to make it easier to premute X(1) etc. See https://github.com/rafaqz/DimensionalData.jl/issues/1144
_plus(x::D,y::D) where D <: Rasters.Dimension{<:Real} = D(x.val + y.val)
function _get_prec_days(force = false)
    climatevars = get_terraclimate_raw((:tavg,))
    climatevars2 = set(set(climatevars, Ti => Intervals(Start())), Ti => Regular(Month(1)))

    filepath = "data/monthly_prec_days.nc"

    if !isfile(filepath) || force
        prec_filepath = _maybe_download_eobs()
        r = Raster(prec_filepath; lazy = true)
        daily_prec = read(crop(r; to = climatevars2))[Ti = 1:End()-1]

        monthly_prec_days = Rasters.combine(RasterSeries((
            dropdims.(count.(
            x -> ismissing(x) ? false : x > 0.2, Rasters.groupby(daily_prec, Ti => yearmonth), 
            dims = :Ti);dims = :Ti)
        )))
        monthly_prec_days = mask(monthly_prec_days; with = daily_prec[Ti = 1])

        monthly_prec_days_res = resample(reorder(monthly_prec_days, dims(climatevars)); to = climatevars, mode = :median)
        #### Fix the rasters - some cells are missing so just fill those with the nearest cell

        # Cycle through and make sure no land values are missing.
        for i in 1:2
            # expand by one in all directions to avoid missing values
            for idx in DimIndices(dims(monthly_prec_days_res, (X,Y)))
                if ismissing(monthly_prec_days_res[idx,Ti = 1])
                    for dir in ((X(-1), Y(0)), (X(1), Y(0)), (X(0), Y(-1)), (X(0), Y(1)))
                        neighbor_idx = _plus.(idx, dir)
                        if Base.checkbounds(Bool, monthly_prec_days_res, CartesianIndex(neighbor_idx[1].val, neighbor_idx[2].val, 1))
                            val = monthly_prec_days_res[neighbor_idx, Ti = 1]
                            if !ismissing(val)
                                monthly_prec_days_res[idx] .= monthly_prec_days_res[neighbor_idx]
                                break
                            end
                        end
                    end
                end
            end
        end
        # replace year month tuple with actual date
        monthly_prec_days_res_to_write = set(monthly_prec_days_res, Ti => 
            broadcast(x -> DateTime(x...), lookup(monthly_prec_days_res, :Ti)))
        write(filepath, monthly_prec_days_res_to_write, force = true)
    end
    Raster(filepath)
end
