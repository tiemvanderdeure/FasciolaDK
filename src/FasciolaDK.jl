module FasciolaDK

using CSV, DataFrames, Dates, Statistics, StatsBase
using Rasters, Rasters.Lookups, NCDatasets, GADM

const datadir = "data/slagtedata"#"/home/tvd/K/slagtedata"

export load_registry, datadir
export get_terraclimate

include("registry.jl")
include("climate.jl")
include("utils.jl")

end # module FasciolaDK
