using GibbsSeaWater
using Test

as = 35.0; ct = 10.0; pr = 10.0
println(GibbsSeaWater.gsw_rho(as, ct, pr))

# load test data
include("gsw_check_data.jl")

# test if value is equal to reference upto a precession of the
# variable reference_ca
# non-finite values in reference are ignored
macro gswtest(value,reference)
    return quote
        mask = isfinite.($reference)
        @test $value[mask] ≈ $reference[mask]  atol=$(Symbol(reference,:_ca))
    end
end


# first dimension is the depth and the
# second dimension is the station
# transpose lat_cast and long_cast so that the
# last dimension is the station
lat = lat_cast'
lon = long_cast'

# check gsw functions

c = GibbsSeaWater.gsw_c_from_sp.(sp,t,p)
@testset "Practical Salinity, PSS-78" begin
    @gswtest GibbsSeaWater.gsw_c_from_sp.(sp,t,p) c_from_sp
    @gswtest GibbsSeaWater.gsw_sp_from_c.(c,t,p) sp_from_c
    @gswtest GibbsSeaWater.gsw_sp_from_sk.(sk) sp_from_sk
end

@testset "Absolute Salinity, Preformed Salinity and Conservative Temperature" begin
    @gswtest GibbsSeaWater.gsw_sa_from_sp.(sp,p,lon,lat) sa_from_sp
    @gswtest GibbsSeaWater.gsw_sstar_from_sp.(sp,p,lon,lat) sstar_from_sp
    @gswtest GibbsSeaWater.gsw_ct_from_t.(sa,t,p) ct_from_t
end

sr = GibbsSeaWater.gsw_sr_from_sp.(sp)
#issue see issue #3
#z = GibbsSeaWater.gsw_z_from_p.(p,lat,0.,0.)
sstar = GibbsSeaWater.gsw_sstar_from_sa.(sa,p,lon,lat)
pt = GibbsSeaWater.gsw_pt_from_ct.(sa,ct)
entropy = GibbsSeaWater.gsw_entropy_from_pt.(sa,pt)
@testset "Other conversions between Temperatures, Salinities, Entropy, Pressure and Height" begin
	@gswtest GibbsSeaWater.gsw_deltasa_from_sp.(sp,p,lon,lat)  deltasa_from_sp
	@gswtest GibbsSeaWater.gsw_sr_from_sp.(sp) sr_from_sp
	@gswtest GibbsSeaWater.gsw_sp_from_sr.(sr) sp_from_sr
	@gswtest GibbsSeaWater.gsw_sp_from_sa.(sa,p,lon,lat) sp_from_sa
	@gswtest GibbsSeaWater.gsw_sstar_from_sa.(sa,p,lon,lat) sstar_from_sa
	@gswtest GibbsSeaWater.gsw_sa_from_sstar.(sstar,p,lon,lat)  sa_from_sstar
	@gswtest GibbsSeaWater.gsw_sp_from_sstar.(sstar,p,lon,lat)  sp_from_sstar
	@gswtest GibbsSeaWater.gsw_pt_from_ct.(sa,ct) pt_from_ct
	@gswtest GibbsSeaWater.gsw_t_from_ct.(sa,ct,p) t_from_ct
	@gswtest GibbsSeaWater.gsw_ct_from_pt.(sa,pt) ct_from_pt
	@gswtest GibbsSeaWater.gsw_pt0_from_t.(sa,t,p) pt0_from_t
	@gswtest GibbsSeaWater.gsw_pt_from_t.(sa,t,p,pref[1]) pt_from_t
#	@gswtest GibbsSeaWater.gsw_z_from_p.(p,lat,0.,0.) z_from_p
#	@gswtest GibbsSeaWater.gsw_p_from_z.(z,lat,0.,0.) p_from_z
	@gswtest GibbsSeaWater.gsw_entropy_from_pt.(sa,pt) entropy_from_pt
	@gswtest GibbsSeaWater.gsw_pt_from_entropy.(sa,entropy) pt_from_entropy
	@gswtest GibbsSeaWater.gsw_ct_from_entropy.(sa,entropy) ct_from_entropy
	@gswtest GibbsSeaWater.gsw_entropy_from_ct.(sa,ct) entropy_from_ct
	@gswtest GibbsSeaWater.gsw_entropy_from_t.(sa,t,p) entropy_from_t
	@gswtest GibbsSeaWater.gsw_adiabatic_lapse_rate_from_ct.(sa,ct,p)  adiabatic_lapse_rate_from_ct
end
