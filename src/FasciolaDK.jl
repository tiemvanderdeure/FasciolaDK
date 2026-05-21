module FasciolaDK

using CSV, DataFrames, Dates, Statistics, StatsBase
using Rasters, Rasters.Lookups, NCDatasets, GADM
import Printf: @sprintf

const datadir = "/home/tvd/K/FasciolaDK/slagtedata"

export load_registry, datadir
export get_terraclimate

include("registry.jl")
include("climate.jl")
include("plots.jl")

end # module FasciolaDK
