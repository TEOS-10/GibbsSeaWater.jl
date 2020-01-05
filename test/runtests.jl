# adapted from
# https://github.com/TEOS-10/GSW-C/blob/95a9caeb128e4793b538b5ded889c83c913584fb/gsw_check_functions.c
#
# Copyright (c) 2011, SCOR/IAPSO WG127 (Scientific Committee on Oceanic Research/
# International Association for the Physical Sciences of the Oceans, Working Group 127).
# https://github.com/TEOS-10/GSW-C/blob/95a9caeb128e4793b538b5ded889c83c913584fb/LICENSE

using GibbsSeaWater
using Test

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

macro gswtest(value,reference...)
    return Expr(:block, [
    quote
        @test map(vab -> vab[$ii],$value) ≈ $(reference[ii]) nans=true atol=$(Symbol(reference[ii],:_ca))
    end for ii in 1:length(reference)]...)
end


# first dimension is the depth and the
# second dimension is the station
# transpose lat_cast and long_cast so that the
# last dimension is the station
lat = lat_cast'
lon = long_cast'

# check gsw functions
c = gsw_c_from_sp.(sp,t,p)
@testset "Practical Salinity, PSS-78" begin
    @gswtest gsw_c_from_sp.(sp,t,p) c_from_sp
    @gswtest gsw_sp_from_c.(c,t,p) sp_from_c
    @gswtest gsw_sp_from_sk.(sk) sp_from_sk
end

@testset "Absolute Salinity, Preformed Salinity and Conservative Temperature" begin
    @gswtest gsw_sa_from_sp.(sp,p,lon,lat) sa_from_sp
    @gswtest gsw_sstar_from_sp.(sp,p,lon,lat) sstar_from_sp
    @gswtest gsw_ct_from_t.(sa,t,p) ct_from_t
end

sr = gsw_sr_from_sp.(sp)
#issue see issue #3
#z = gsw_z_from_p.(p,lat,0.,0.)
sstar = gsw_sstar_from_sa.(sa,p,lon,lat)
pt = gsw_pt_from_ct.(sa,ct)
entropy = gsw_entropy_from_pt.(sa,pt)
@testset "Other conversions between Temperatures, Salinities, Entropy, Pressure and Height" begin
	@gswtest gsw_deltasa_from_sp.(sp,p,lon,lat)  deltasa_from_sp
	@gswtest gsw_sr_from_sp.(sp) sr_from_sp
	@gswtest gsw_sp_from_sr.(sr) sp_from_sr
	@gswtest gsw_sp_from_sa.(sa,p,lon,lat) sp_from_sa
	@gswtest gsw_sstar_from_sa.(sa,p,lon,lat) sstar_from_sa
	@gswtest gsw_sa_from_sstar.(sstar,p,lon,lat)  sa_from_sstar
	@gswtest gsw_sp_from_sstar.(sstar,p,lon,lat)  sp_from_sstar
	@gswtest gsw_pt_from_ct.(sa,ct) pt_from_ct
	@gswtest gsw_t_from_ct.(sa,ct,p) t_from_ct
	@gswtest gsw_ct_from_pt.(sa,pt) ct_from_pt
	@gswtest gsw_pt0_from_t.(sa,t,p) pt0_from_t
	@gswtest gsw_pt_from_t.(sa,t,p,pref[1]) pt_from_t
#	@gswtest gsw_z_from_p.(p,lat,0.,0.) z_from_p
#	@gswtest gsw_p_from_z.(z,lat,0.,0.) p_from_z
	@gswtest gsw_entropy_from_pt.(sa,pt) entropy_from_pt
	@gswtest gsw_pt_from_entropy.(sa,entropy) pt_from_entropy
	@gswtest gsw_ct_from_entropy.(sa,entropy) ct_from_entropy
	@gswtest gsw_entropy_from_ct.(sa,ct) entropy_from_ct
	@gswtest gsw_entropy_from_t.(sa,t,p) entropy_from_t
	@gswtest gsw_adiabatic_lapse_rate_from_ct.(sa,ct,p)  adiabatic_lapse_rate_from_ct
end



@testset "Specific Volume, Density and Enthalpy" begin
    @gswtest gsw_specvol.(sa,ct,p) specvol
	@gswtest gsw_alpha.(sa,ct,p) alpha
	@gswtest gsw_beta.(sa,ct,p) beta
	@gswtest gsw_alpha_on_beta.(sa,ct,p) alpha_on_beta
    @gswtest gsw_specvol_alpha_beta.(sa,ct,p) v_vab alpha_vab beta_vab
	@gswtest gsw_specvol_first_derivatives.(sa,ct,p) v_sa v_ct v_p
	@gswtest gsw_specvol_second_derivatives.(sa,ct,p) v_sa_sa v_sa_ct v_ct_ct v_sa_p v_ct_p
#=
	test_sub2(specvol_first_derivatives_wrt_enthalpy,
	    (sa[i],ct[i],p[i],VALS2),val1,v_sa_wrt_h, val2,v_h);
	test_sub3(specvol_second_derivatives_wrt_enthalpy,
	    (sa[i],ct[i],p[i],VALS3),
	    val1,v_sa_sa_wrt_h, val2,v_sa_h, val3,v_h_h);
	@gswtest gsw_specvol_anom_standard.(sa[i],ct[i],p[i]), value,
	    specvol_anom_standard);
	@gswtest gsw_rho.(sa[i],ct[i],p[i]), rho,rho);
	test_sub3(rho_alpha_beta, (sa[i],ct[i],p[i],VALS3),
	    val1,rho_rab,val2,alpha_rab,val3,beta_rab);
	test_sub3(rho_first_derivatives,(sa[i],ct[i],p[i],VALS3),
	    val1,rho_sa,val2,rho_ct,val3,rho_p);
	test_sub5(rho_second_derivatives,
	    (sa[i],ct[i],p[i],VALS5),val1,rho_sa_sa,val2,rho_sa_ct,
	    val3,rho_ct_ct,val4,rho_sa_p,val5, rho_ct_p);
	test_sub2(rho_first_derivatives_wrt_enthalpy,
	    (sa[i],ct[i],p[i],VALS2),val1,rho_sa_wrt_h,val2,rho_h);
	test_sub3(rho_second_derivatives_wrt_enthalpy,
		(sa[i],ct[i],p[i],VALS3),val1,rho_sa_sa_wrt_h,
		val2,rho_sa_h, val3,rho_h_h);
	@gswtest gsw_sigma0.(sa[i],ct[i]), value,sigma0);
	@gswtest gsw_sigma1.(sa[i],ct[i]), value,sigma1);
	@gswtest gsw_sigma2.(sa[i],ct[i]), value,sigma2);
	@gswtest gsw_sigma3.(sa[i],ct[i]), value,sigma3);
	@gswtest gsw_sigma4.(sa[i],ct[i]), value,sigma4);
	@gswtest gsw_sound_speed.(sa[i],ct[i],p[i]), value,sound_speed);
	@gswtest gsw_kappa.(sa[i],ct[i],p[i]), value,kappa);
	@gswtest gsw_cabbeling.(sa[i],ct[i],p[i]), value,cabbeling);
	@gswtest gsw_thermobaric.(sa[i],ct[i],p[i]), value,thermobaric);
	@gswtest gsw_sa_from_rho.(rho[i],ct[i],p[i]), value,sa_from_rho);
	test_sub1(ct_from_rho, (rho[i],sa[i],p[i],VALS1,NULL),val1,ct_from_rho);
	@gswtest gsw_ct_maxdensity.(sa[i],p[i]), value,ct_maxdensity);
	@gswtest gsw_internal_energy.(sa[i],ct[i],p[i]), value,internal_energy);
	@gswtest gsw_enthalpy.(sa[i],ct[i],p[i]), h,enthalpy);
	@gswtest gsw_enthalpy_diff,
	    (sa[i],ct[i],p_shallow[i],p_deep[i]), value,enthalpy_diff);
	@gswtest gsw_ct_from_enthalpy.(sa[i],h[i],p[i]), value,ct_from_enthalpy);
	@gswtest gsw_dynamic_enthalpy.(sa[i],ct[i],p[i]), value,dynamic_enthalpy);
	test_sub2(enthalpy_first_derivatives, (sa[i],ct[i],p[i],VALS2),
	    val1,h_sa, val2,h_ct);
	test_sub3(enthalpy_second_derivatives,(sa[i],ct[i],p[i],VALS3),
	    val1,h_sa_sa,val2,h_sa_ct, val3,h_ct_ct);
=#
end
#=
	@testset "Derivatives of entropy, CT and pt" begin

	test_sub2(ct_first_derivatives, (sa[i],pt[i],VALS2),
	    val1,ct_sa, val2,ct_pt);
	test_sub3(ct_second_derivatives, (sa[i],pt[i],VALS3),
	    val1,ct_sa_sa, val2,ct_sa_pt, val3,ct_pt_pt);
	test_sub2(entropy_first_derivatives, (sa[i],ct[i],VALS2),
	    val1,eta_sa, val2,eta_ct);
	test_sub3(entropy_second_derivatives, (sa[i],ct[i],VALS3),
	    val1,eta_sa_sa, val2,eta_sa_ct, val3,eta_ct_ct);
	test_sub2(pt_first_derivatives, (sa[i],ct[i],VALS2),
	    val1,pt_sa, val2,pt_ct);
	test_sub3(pt_second_derivatives, (sa[i],ct[i],VALS3),
	    val1,pt_sa_sa,val2,pt_sa_ct,val3,pt_ct_ct);
end

	@testset "Freezing temperatures" begin

	saturation_fraction = 0.5;

	@gswtest gsw_ct_freezing,(sa[i],p[i],saturation_fraction),
	    ctf,ct_freezing);
	@gswtest gsw_ct_freezing_poly.(sa[i],p[i],saturation_fraction), ctf_poly,
	    ct_freezing_poly);
	@gswtest gsw_t_freezing.(sa[i],p[i],saturation_fraction), tf,
	    t_freezing);
	@gswtest gsw_t_freezing_poly.(sa[i],p[i],saturation_fraction), tf_poly,
	    t_freezing_poly);
	@gswtest gsw_pot_enthalpy_ice_freezing.(sa[i],p[i]), value,
	    pot_enthalpy_ice_freezing);
	@gswtest gsw_pot_enthalpy_ice_freezing_poly.(sa[i],p[i]), value,
	    pot_enthalpy_ice_freezing_poly);
	@gswtest gsw_sa_freezing_from_ct.(ctf[i],p[i],saturation_fraction),value,
	    sa_freezing_from_ct);
	@gswtest gsw_sa_freezing_from_ct_poly,
	    (ctf_poly[i],p[i],saturation_fraction),value,
	    sa_freezing_from_ct_poly);
	@gswtest gsw_sa_freezing_from_t.(tf[i],p[i],saturation_fraction),value,
	    sa_freezing_from_t);
	@gswtest gsw_sa_freezing_from_t_poly,
	    (tf_poly[i],p[i],saturation_fraction),value,
	    sa_freezing_from_t_poly);
	test_sub2(ct_freezing_first_derivatives,
	    (sa[i],p[i],saturation_fraction,VALS2),
		val1,ctfreezing_sa,val2,ctfreezing_p);
	test_sub2(ct_freezing_first_derivatives_poly,
	    (sa[i],p[i],saturation_fraction,VALS2),val1,ctfreezing_sa_poly,
		val2,ctfreezing_p_poly);
	test_sub2(t_freezing_first_derivatives,
	    (sa[i],p[i],saturation_fraction,VALS2),
	    val1,tfreezing_sa,val2,tfreezing_p);
	test_sub2(t_freezing_first_derivatives_poly,
	    (sa[i],p[i],saturation_fraction,VALS2),
	    val1,tfreezing_sa_poly,val2,tfreezing_p_poly);
	test_sub2(pot_enthalpy_ice_freezing_first_derivatives,
		(sa[i],p[i],VALS2),
		val1, pot_enthalpy_ice_freezing_sa,
		val2, pot_enthalpy_ice_freezing_p);
	test_sub2(pot_enthalpy_ice_freezing_first_derivatives_poly,
		(sa[i],p[i],VALS2),
		val1, pot_enthalpy_ice_freezing_sa_poly,
		val2, pot_enthalpy_ice_freezing_p_poly);

	end

	@testset "Isobaric Melting Enthalpy and Isobaric Evaporation Enthalpy" begin

	@gswtest gsw_latentheat_melting.(sa[i],p[i]), value,latentheat_melting);
	@gswtest gsw_latentheat_evap_ct.(sa[i],ct[i]), value,latentheat_evap_ct);
	@gswtest gsw_latentheat_evap_t.(sa[i],t[i]), value,latentheat_evap_t);

	end

	@testset "Planet Earth properties" begin

	@gswtest gsw_grav.(lat[i],p[i]), value,grav);

	end

	@testset "Density and enthalpy in terms of CT, derived from the exact Gibbs function" begin

	@gswtest gsw_enthalpy_ct_exact,(sa[i],ct[i],p[i]),value,enthalpy_ct_exact);
	test_sub2(enthalpy_first_derivatives_ct_exact, (sa[i],ct[i],p[i],VALS2),
	    val1,h_sa_ct_exact,val2,h_ct_ct_exact);
	test_sub3(enthalpy_second_derivatives_ct_exact, (sa[i],ct[i],p[i],VALS3),
	    val1,h_sa_sa_ct_exact, val2,h_sa_ct_ct_exact,val3,h_ct_ct_ct_exact);

	end

	@testset "Basic thermodynamic properties in terms of in-situ t, based on the exact Gibbs function" begin

	@gswtest gsw_rho_t_exact.(sa[i],t[i],p[i]),value,rho_t_exact);
	@gswtest gsw_pot_rho_t_exact.(sa[i],t[i],p[i],pref[0]),value,
	    pot_rho_t_exact);
	@gswtest gsw_alpha_wrt_t_exact.(sa[i],t[i],p[i]),value,alpha_wrt_t_exact);
	@gswtest gsw_beta_const_t_exact.(sa[i],t[i],p[i]),value,
	    beta_const_t_exact);
	@gswtest gsw_specvol_t_exact.(sa[i],t[i],p[i]),value,specvol_t_exact);
	@gswtest gsw_sound_speed_t_exact.(sa[i],t[i],p[i]),value,
	    sound_speed_t_exact);
	@gswtest gsw_kappa_t_exact.(sa[i],t[i],p[i]),value,kappa_t_exact);
	@gswtest gsw_enthalpy_t_exact.(sa[i],t[i],p[i]),value,enthalpy_t_exact);
	test_sub3(ct_first_derivatives_wrt_t_exact, (sa[i],t[i],p[i],VALS3),
	    val1,ct_sa_wrt_t, val2,ct_t_wrt_t,val3,ct_p_wrt_t);
	@gswtest gsw_chem_potential_water_t_exact.(sa[i],t[i],p[i]),value,
	    chem_potential_water_t_exact);
	@gswtest gsw_t_deriv_chem_potential_water_t_exact.(sa[i],t[i],p[i]),
	    value,t_deriv_chem_potential_water_t_exact);
	@gswtest gsw_dilution_coefficient_t_exact.(sa[i],t[i],p[i]),value,
	    dilution_coefficient_t_exact);

	end

	@testset "Library functions of the GSW Toolbox" begin

	@gswtest gsw_deltasa_atlas.(p[i],lon[i],lat[i]),value,deltasa_atlas);
	@gswtest gsw_fdelta.(p[i],lon[i],lat[i]),value,fdelta);

	end

	@testset "Water column properties, based on the 75-term polynomial " *
	    "for specific volume" begin

	count	= cast_mpres_m*cast_mpres_n;
	for (j = 0; j<cast_mpres_n; j++) {
	    k = j*cast_m; l = j*cast_mpres_m;
    	    gsw_nsquared(&sa[k],&ct[k],&p[k],&lat[k],cast_m,&val1[l],&val2[l]);
	}
	check_accuracy("nsquared",n2_ca,"n2",count, val1, n2);
	check_accuracy("nsquared",p_mid_n2_ca,"p_mid_n2",count, val2, p_mid_n2);

	for (j = 0; j<cast_mpres_n; j++) {
	    k = j*cast_m; l = j*cast_mpres_m;
	    gsw_turner_rsubrho(&sa[k],&ct[k],&p[k],cast_m,&val1[l],&val2[l],
				&val3[l]);
	}
	check_accuracy("turner_rsubrho",tu_ca,"tu",count, val1, tu);
	check_accuracy("rsubrhorner_rsubrho",rsubrho_ca,"rsubrho",count, val2,
		rsubrho);
	check_accuracy("p_mid_tursrrner_rsubrho",p_mid_tursr_ca,"p_mid_tursr",
		count, val3, p_mid_tursr);

	for (j = 0; j<cast_mpres_n; j++) {
	    k = j*cast_m; l = j*cast_mpres_m;
	    gsw_ipv_vs_fnsquared_ratio(&sa[k],&ct[k],&p[k],pref[0],cast_m,
		&val1[l], &val2[l]);
	}
	check_accuracy("ipv_vs_fnsquared_ratio",ipvfn2_ca,"ipvfn2",count,
		val1, ipvfn2);
	check_accuracy("ipv_vs_fnsquared_ratio",p_mid_ipvfn2_ca,"p_mid_ipvfn2",
		count, val2, p_mid_ipvfn2);

	for (j = 0; j<cast_mpres_n; j++) {
	    k = j*cast_m;
	    for (n=0; n<cast_m; n++)
		if (isnan(sa[k+n]) || fabs(sa[k+n]) >= GSW_ERROR_LIMIT)
		    break;
	    if (gsw_geo_strf_dyn_height(&sa[k],&ct[k],&p[k],pref[0],n,
		&val1[k]) == NULL)
		printf("geo_strf_dyn_height returned NULL.\n");
	}
	check_accuracy("geo_strf_dyn_height",geo_strf_dyn_height_ca,
		"geo_strf_dyn_height",count, val1, geo_strf_dyn_height);

	for (j = 0; j<cast_mpres_n; j++) {
	    k = j*cast_m;
	    for (n=0; n<cast_m; n++)
		if (isnan(sa[k+n]) || fabs(sa[k+n]) >= GSW_ERROR_LIMIT)
		    break;
	    gsw_geo_strf_dyn_height_pc(&sa[k],&ct[k],&delta_p[k],n,
		&val1[k], &val2[k]);
	}
	check_accuracy("geo_strf_dyn_height_pc",geo_strf_dyn_height_pc_ca,
		"geo_strf_dyn_height_pc",count, val1, geo_strf_dyn_height_pc);
	check_accuracy("geo_strf_dyn_height_pc",geo_strf_dyn_height_pc_p_mid_ca,
		"geo_strf_dyn_height_pc_p_mid",count, val2,
		geo_strf_dyn_height_pc_p_mid);

	end

	@testset "Thermodynamic properties of ice Ih" begin

	count = cast_ice_m*cast_ice_n;

	@gswtest gsw_rho_ice.(t_seaice[i],p_arctic[i]),value,rho_ice);
	@gswtest gsw_alpha_wrt_t_ice.(t_seaice[i],p_arctic[i]),value,
	    alpha_wrt_t_ice);
	@gswtest gsw_specvol_ice.(t_seaice[i],p_arctic[i]),value,specvol_ice);
	@gswtest gsw_pressure_coefficient_ice.(t_seaice[i],p_arctic[i]),value,
	    pressure_coefficient_ice);
	@gswtest gsw_sound_speed_ice.(t_seaice[i],p_arctic[i]),value,
	    sound_speed_ice);
	@gswtest gsw_kappa_ice.(t_seaice[i],p_arctic[i]),value,kappa_ice);
	@gswtest gsw_kappa_const_t_ice.(t_seaice[i],p_arctic[i]),value,
	    kappa_const_t_ice);
	@gswtest gsw_internal_energy_ice.(t_seaice[i],p_arctic[i]),value,
	    internal_energy_ice);
	@gswtest gsw_enthalpy_ice.(t_seaice[i],p_arctic[i]),value,enthalpy_ice);
	@gswtest gsw_entropy_ice.(t_seaice[i],p_arctic[i]),value,entropy_ice);
	@gswtest gsw_cp_ice.(t_seaice[i],p_arctic[i]),value,cp_ice);
	@gswtest gsw_chem_potential_water_ice,
	    (t_seaice[i],p_arctic[i]),value,chem_potential_water_ice);
	@gswtest gsw_helmholtz_energy_ice.(t_seaice[i],p_arctic[i]),value,
	    helmholtz_energy_ice);
	@gswtest gsw_adiabatic_lapse_rate_ice,
	    (t_seaice[i],p_arctic[i]),value,adiabatic_lapse_rate_ice);
	@gswtest gsw_pt0_from_t_ice,(t_seaice[i],p_arctic[i]),pt0, pt0_from_t_ice);
	@gswtest gsw_pt_from_t_ice.(t_seaice[i],p_arctic[i],pref[0]),value,
	    pt_from_t_ice);
	@gswtest gsw_t_from_pt0_ice.(pt0[i],p_arctic[i]),value,t_from_pt0_ice);
	@gswtest gsw_pot_enthalpy_from_pt_ice.(pt0[i]), h,
	    pot_enthalpy_from_pt_ice);
	@gswtest gsw_pt_from_pot_enthalpy_ice.(h[i]), value,
	    pt_from_pot_enthalpy_ice);
	@gswtest gsw_pot_enthalpy_from_pt_ice_poly.(pt0[i]), h,
	    pot_enthalpy_from_pt_ice_poly);
	@gswtest gsw_pt_from_pot_enthalpy_ice_poly.(h[i]),value,
	    pt_from_pot_enthalpy_ice_poly);

	saturation_fraction = 0.5;

	@gswtest gsw_pressure_freezing_ct,
	    (sa_arctic[i],ct_arctic[i]-1.0,saturation_fraction),value,
	    pressure_freezing_ct);

	end

	@testset "Thermodynamic interaction between ice and seawater" begin

	@gswtest gsw_melting_ice_sa_ct_ratio,
	    (sa_arctic[i],ct_arctic[i],p_arctic[i],t_ice[i]),value,
	    melting_ice_sa_ct_ratio);
	@gswtest gsw_melting_ice_sa_ct_ratio_poly,
	    (sa_arctic[i],ct_arctic[i],p_arctic[i],t_ice[i]),value,
	    melting_ice_sa_ct_ratio_poly);
	@gswtest gsw_melting_ice_equilibrium_sa_ct_ratio,
	    (sa_arctic[i],p_arctic[i]),value,
	    melting_ice_equilibrium_sa_ct_ratio);
	@gswtest gsw_melting_ice_equilibrium_sa_ct_ratio_poly,
	    (sa_arctic[i],p_arctic[i]),value,
	    melting_ice_equilibrium_sa_ct_ratio_poly);
	test_sub2(melting_ice_into_seawater,
	    (sa_arctic[i],ct_arctic[i]+0.1,p_arctic[i],w_ice[i],t_ice[i],VALS3),
	    val1, melting_ice_into_seawater_sa_final,
	    val2, melting_ice_into_seawater_ct_final);
	    /*val3, melting_ice_into_seawater_w_ih);*/
	test_sub3(ice_fraction_to_freeze_seawater,
	    (sa_arctic[i],ct_arctic[i],p_arctic[i],t_ice[i],VALS3),
	    val1, ice_fraction_to_freeze_seawater_sa_freeze,
	    val2, ice_fraction_to_freeze_seawater_ct_freeze,
	    val3, ice_fraction_to_freeze_seawater_w_ih);
	test_sub3(frazil_ratios_adiabatic,
	    (sa_arctic[i],p_arctic[i],w_ice[i],VALS3),
	    val1,dsa_dct_frazil, val2,dsa_dp_frazil,
	    val3,dct_dp_frazil);
	test_sub3(frazil_ratios_adiabatic_poly,
	    (sa_arctic[i],p_arctic[i],w_ice[i],VALS3),
	    val1,dsa_dct_frazil_poly, val2,dsa_dp_frazil_poly,
	    val3,dct_dp_frazil_poly);
	test_sub3(frazil_properties_potential,
	    (sa_bulk[i],h_pot_bulk[i],p_arctic[i], VALS3),
	    val1, frazil_properties_potential_sa_final,
	    val2, frazil_properties_potential_ct_final,
	    val3, frazil_properties_potential_w_ih_final);
	test_sub3(frazil_properties_potential_poly,
	    (sa_bulk[i],h_pot_bulk[i], p_arctic[i],VALS3),
	    val1, frazil_properties_potential_poly_sa_final,
	    val2, frazil_properties_potential_poly_ct_final,
	    val3, frazil_properties_potential_poly_w_ih_final);
	test_sub3(frazil_properties,
	    (sa_bulk[i],h_bulk[i],p_arctic[i],VALS3),
	    val1,frazil_properties_sa_final,
	    val2,frazil_properties_ct_final,
	    val3,frazil_properties_w_ih_final);

	end

	@testset "Thermodynamic interaction between seaice and seawater" begin

	@gswtest gsw_melting_seaice_sa_ct_ratio,
	    (sa_arctic[i],ct_arctic[i],p_arctic[i], sa_seaice[i],t_seaice[i]),
	    value,melting_seaice_sa_ct_ratio);
	@gswtest gsw_melting_seaice_sa_ct_ratio_poly,
	    (sa_arctic[i],ct_arctic[i],p_arctic[i], sa_seaice[i],t_seaice[i]),
	    value,melting_seaice_sa_ct_ratio_poly);
	@gswtest gsw_melting_seaice_equilibrium_sa_ct_ratio,
	    (sa_arctic[i],p_arctic[i]),value,
	    melting_seaice_equilibrium_sa_ct_ratio);
	@gswtest gsw_melting_seaice_equilibrium_sa_ct_ratio_poly,
	    (sa_arctic[i],p_arctic[i]),value,
	    melting_seaice_equilibrium_sa_ct_ratio_poly);
	test_sub2(melting_seaice_into_seawater, (sa_arctic[i],ct_arctic[i],
	    p_arctic[i], w_seaice[i],sa_seaice[i],t_seaice[i],VALS2),
	    val1, melting_seaice_into_seawater_sa_final,
	    val2, melting_seaice_into_seawater_ct_final);
	test_sub3(seaice_fraction_to_freeze_seawater,(sa_arctic[i],ct_arctic[i],
	    p_arctic[i], sa_seaice[i],t_seaice[i],VALS3),
	    val1, seaice_fraction_to_freeze_seawater_sa_freeze,
	    val2, seaice_fraction_to_freeze_seawater_ct_freeze,
	    val3, seaice_fraction_to_freeze_seawater_w_ih);

end
=#
