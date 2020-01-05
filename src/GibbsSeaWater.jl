module GibbsSeaWater

deps = joinpath(dirname(dirname(pathof(GibbsSeaWater))), "deps", "deps.jl")
isfile(deps) ? include(deps) : error("GibbsSeaWater is not properly installed")

include("gen_gswteos_h.jl")
include("gen_gswteos10.jl")

function gsw_specvol_alpha_beta(sa::Float64, ct::Float64, p::Float64)
    specvol = Ref{Float64}()
    alpha = Ref{Float64}()
    beta = Ref{Float64}()
    gsw_specvol_alpha_beta(sa,ct,p,specvol,alpha,beta)
    return specvol[],alpha[],beta[]
end

function gsw_specvol_first_derivatives(sa::Cdouble, ct::Cdouble, p::Cdouble)
    v_sa = Ref{Cdouble}()
    v_ct = Ref{Cdouble}()
    v_p = Ref{Cdouble}()
    gsw_specvol_first_derivatives(sa, ct, p, v_sa, v_ct, v_p)
    return v_sa[], v_ct[], v_p[]
end

function gsw_specvol_second_derivatives(sa::Cdouble, ct::Cdouble, p::Cdouble)
    v_sa_sa = Ref{Cdouble}()
    v_sa_ct = Ref{Cdouble}()
    v_ct_ct = Ref{Cdouble}()
    v_sa_p = Ref{Cdouble}()
    v_ct_p = Ref{Cdouble}()

    gsw_specvol_second_derivatives(sa, ct, p, v_sa_sa, v_sa_ct, v_ct_ct, v_sa_p, v_ct_p)
    return v_sa_sa[], v_sa_ct[], v_ct_ct[], v_sa_p[], v_ct_p[]
end

export gsw_add_barrier
export gsw_add_mean
export gsw_adiabatic_lapse_rate_from_ct
export gsw_adiabatic_lapse_rate_ice
export gsw_alpha
export gsw_alpha_on_beta
export gsw_alpha_wrt_t_exact
export gsw_alpha_wrt_t_ice
export gsw_beta_const_t_exact
export gsw_beta
export gsw_cabbeling
export gsw_c_from_sp
export gsw_chem_potential_water_ice
export gsw_chem_potential_water_t_exact
export gsw_cp_ice
export gsw_cp_t_exact
export gsw_ct_first_derivatives
export gsw_ct_first_derivatives_wrt_t_exact
export gsw_ct_freezing
export gsw_ct_freezing_first_derivatives
export gsw_ct_freezing_first_derivatives_poly
export gsw_ct_freezing_poly
export gsw_ct_from_enthalpy
export gsw_ct_from_enthalpy_exact
export gsw_ct_from_entropy
export gsw_ct_from_pt
export gsw_ct_from_rho
export gsw_ct_from_t
export gsw_ct_maxdensity
export gsw_ct_second_derivatives
export gsw_deltasa_atlas
export gsw_deltasa_from_sp
export gsw_dilution_coefficient_t_exact
export gsw_dynamic_enthalpy
export gsw_enthalpy_ct_exact
export gsw_enthalpy_diff
export gsw_enthalpy
export gsw_enthalpy_first_derivatives_ct_exact
export gsw_enthalpy_first_derivatives
export gsw_enthalpy_ice
export gsw_enthalpy_second_derivatives_ct_exact
export gsw_enthalpy_second_derivatives
export gsw_enthalpy_sso_0
export gsw_enthalpy_t_exact
export gsw_entropy_first_derivatives
export gsw_entropy_from_ct
export gsw_entropy_from_pt
export gsw_entropy_from_t
export gsw_entropy_ice
export gsw_entropy_part
export gsw_entropy_part_zerop
export gsw_entropy_second_derivatives
export gsw_fdelta
export gsw_frazil_properties
export gsw_frazil_properties_potential
export gsw_frazil_properties_potential_poly
export gsw_frazil_ratios_adiabatic
export gsw_frazil_ratios_adiabatic_poly
export gsw_geo_strf_dyn_height
export gsw_geo_strf_dyn_height_1
export gsw_geo_strf_dyn_height_pc
export gsw_gibbs_ice
export gsw_gibbs_ice_part_t
export gsw_gibbs_ice_pt0
export gsw_gibbs_ice_pt0_pt0
export gsw_gibbs
export gsw_gibbs_pt0_pt0
export gsw_grav
export gsw_helmholtz_energy_ice
export gsw_hill_ratio_at_sp2
export gsw_ice_fraction_to_freeze_seawater
export gsw_internal_energy
export gsw_internal_energy_ice
export gsw_ipv_vs_fnsquared_ratio
export gsw_kappa_const_t_ice
export gsw_kappa
export gsw_kappa_ice
export gsw_kappa_t_exact
export gsw_latentheat_evap_ct
export gsw_latentheat_evap_t
export gsw_latentheat_melting
export gsw_linear_interp_sa_ct
export gsw_melting_ice_equilibrium_sa_ct_ratio
export gsw_melting_ice_equilibrium_sa_ct_ratio_poly
export gsw_melting_ice_into_seawater
export gsw_melting_ice_sa_ct_ratio
export gsw_melting_ice_sa_ct_ratio_poly
export gsw_melting_seaice_equilibrium_sa_ct_ratio
export gsw_melting_seaice_equilibrium_sa_ct_ratio_poly
export gsw_melting_seaice_into_seawater
export gsw_melting_seaice_sa_ct_ratio
export gsw_melting_seaice_sa_ct_ratio_poly
export gsw_nsquared
export gsw_o2sol
export gsw_o2sol_sp_pt
export gsw_pot_enthalpy_from_pt_ice
export gsw_pot_enthalpy_from_pt_ice_poly
export gsw_pot_enthalpy_ice_freezing
export gsw_pot_enthalpy_ice_freezing_first_derivatives
export gsw_pot_enthalpy_ice_freezing_first_derivatives_poly
export gsw_pot_enthalpy_ice_freezing_poly
export gsw_pot_rho_t_exact
export gsw_pressure_coefficient_ice
export gsw_pressure_freezing_ct
export gsw_pt0_cold_ice_poly
export gsw_pt0_from_t
export gsw_pt0_from_t_ice
export gsw_pt_first_derivatives
export gsw_pt_from_ct
export gsw_pt_from_entropy
export gsw_pt_from_pot_enthalpy_ice
export gsw_pt_from_pot_enthalpy_ice_poly_dh
export gsw_pt_from_pot_enthalpy_ice_poly
export gsw_pt_from_t
export gsw_pt_from_t_ice
export gsw_pt_second_derivatives
export gsw_rho_alpha_beta
export gsw_rho
export gsw_rho_first_derivatives
export gsw_rho_first_derivatives_wrt_enthalpy
export gsw_rho_ice
export gsw_rho_second_derivatives
export gsw_rho_second_derivatives_wrt_enthalpy
export gsw_rho_t_exact
export gsw_rr68_interp_sa_ct
export gsw_saar
export gsw_sa_freezing_estimate
export gsw_sa_freezing_from_ct
export gsw_sa_freezing_from_ct_poly
export gsw_sa_freezing_from_t
export gsw_sa_freezing_from_t_poly
export gsw_sa_from_rho
export gsw_sa_from_sp_baltic
export gsw_sa_from_sp
export gsw_sa_from_sstar
export gsw_sa_p_inrange
export gsw_seaice_fraction_to_freeze_seawater
export gsw_sigma0
export gsw_sigma1
export gsw_sigma2
export gsw_sigma3
export gsw_sigma4
export gsw_sound_speed
export gsw_sound_speed_ice
export gsw_sound_speed_t_exact
export gsw_specvol_alpha_beta
export gsw_specvol_anom_standard
export gsw_specvol
export gsw_specvol_first_derivatives
export gsw_specvol_first_derivatives_wrt_enthalpy
export gsw_specvol_ice
export gsw_specvol_second_derivatives
export gsw_specvol_second_derivatives_wrt_enthalpy
export gsw_specvol_sso_0
export gsw_specvol_t_exact
export gsw_sp_from_c
export gsw_sp_from_sa_baltic
export gsw_sp_from_sa
export gsw_sp_from_sk
export gsw_sp_from_sr
export gsw_sp_from_sstar
export gsw_spiciness0
export gsw_spiciness1
export gsw_spiciness2
export gsw_sr_from_sp
export gsw_sstar_from_sa
export gsw_sstar_from_sp
export gsw_t_deriv_chem_potential_water_t_exact
export gsw_t_freezing
export gsw_t_freezing_first_derivatives_poly
export gsw_t_freezing_first_derivatives
export gsw_t_freezing_poly
export gsw_t_from_ct
export gsw_t_from_pt0_ice
export gsw_thermobaric
export gsw_turner_rsubrho
export gsw_util_indx
export gsw_util_interp1q_int
export gsw_util_linear_interp
export gsw_util_sort_real
export gsw_util_xinterp1
export gsw_util_pchip_interp
export gsw_z_from_p
export gsw_p_from_z

end # module
