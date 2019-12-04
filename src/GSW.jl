__precompile__()

module GSW

deps = joinpath(dirname(dirname(pathof(GSW))), "deps", "deps.jl")
isfile(deps) ? include(deps) : error("GSW is not properly installed")

include("gen_gswteos_h.jl")
include("gen_gswteos10.jl")
# package code goes here



end # module
