using GibbsSeaWater
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
    as = 35.0; ct = 10.0; pr = 10.0
    println(GibbsSeaWater.gsw_rho(as, ct, pr))
end

# load test data
include("gsw_check_data.jl")

# check gsw functions

@testset "Practical Salinity, PSS-78" begin
    mask = isfinite.(c_from_sp)
    @test GibbsSeaWater.gsw_c_from_sp.(sp,t,p)[mask] ≈ c_from_sp[mask] atol=c_from_sp_ca

end
