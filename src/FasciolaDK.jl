module FasciolaDK

using CSV, DataFrames, Dates, Statistics, StatsBase
using Rasters, Rasters.Lookups, NCDatasets, GADM

const datadir = "/home/tvd/K/slagtedata"

export load_registry, datadir
export get_terraclimate
export vsup_colormatrix, val_u_to_color, vsup_legend

include("registry.jl")
include("climate.jl")
include("utils.jl")
include("vsup.jl")
include("plots.jl")

end # module FasciolaDK
