__precompile__()

module GibbsSeaWater

deps = joinpath(dirname(dirname(pathof(GibbsSeaWater))), "deps", "deps.jl")
isfile(deps) ? include(deps) : error("GibbsSeaWater is not properly installed")

include("gen_gswteos_h.jl")
include("gen_gswteos10.jl")
# package code goes here



end # module
