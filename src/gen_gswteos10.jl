# Julia wrapper for header: src/gswteos-10.h
# Automatically generated using Clang.jl wrap_c, version 0.0.0


function gsw_add_barrier(input_data, lon::Cdouble, lat::Cdouble, long_grid::Cdouble, lat_grid::Cdouble, dlong_grid::Cdouble, dlat_grid::Cdouble, output_data)
    ccall((:gsw_add_barrier, libgswteos), Nothing, (Ptr{Cdouble}, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}), input_data, lon, lat, long_grid, lat_grid, dlong_grid, dlat_grid, output_data)
end

function gsw_add_mean(data_in, data_out)
    ccall((:gsw_add_mean, libgswteos), Nothing, (Ptr{Cdouble}, Ptr{Cdouble}), data_in, data_out)
end

function gsw_adiabatic_lapse_rate_from_ct(sa::Cdouble, ct::Cdouble, p::Cdouble)
    ccall((:gsw_adiabatic_lapse_rate_from_ct, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, ct, p)
end

function gsw_adiabatic_lapse_rate_ice(t::Cdouble, p::Cdouble)
    ccall((:gsw_adiabatic_lapse_rate_ice, libgswteos), Cdouble, (Cdouble, Cdouble), t, p)
end

function gsw_alpha(sa::Cdouble, ct::Cdouble, p::Cdouble)
    ccall((:gsw_alpha, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, ct, p)
end

function gsw_alpha_on_beta(sa::Cdouble, ct::Cdouble, p::Cdouble)
    ccall((:gsw_alpha_on_beta, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, ct, p)
end

function gsw_alpha_wrt_t_exact(sa::Cdouble, t::Cdouble, p::Cdouble)
    ccall((:gsw_alpha_wrt_t_exact, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, t, p)
end

function gsw_alpha_wrt_t_ice(t::Cdouble, p::Cdouble)
    ccall((:gsw_alpha_wrt_t_ice, libgswteos), Cdouble, (Cdouble, Cdouble), t, p)
end

function gsw_beta_const_t_exact(sa::Cdouble, t::Cdouble, p::Cdouble)
    ccall((:gsw_beta_const_t_exact, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, t, p)
end

function gsw_beta(sa::Cdouble, ct::Cdouble, p::Cdouble)
    ccall((:gsw_beta, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, ct, p)
end

function gsw_cabbeling(sa::Cdouble, ct::Cdouble, p::Cdouble)
    ccall((:gsw_cabbeling, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, ct, p)
end

function gsw_c_from_sp(sp::Cdouble, t::Cdouble, p::Cdouble)
    ccall((:gsw_c_from_sp, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sp, t, p)
end

function gsw_chem_potential_water_ice(t::Cdouble, p::Cdouble)
    ccall((:gsw_chem_potential_water_ice, libgswteos), Cdouble, (Cdouble, Cdouble), t, p)
end

function gsw_chem_potential_water_t_exact(sa::Cdouble, t::Cdouble, p::Cdouble)
    ccall((:gsw_chem_potential_water_t_exact, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, t, p)
end

function gsw_cp_ice(t::Cdouble, p::Cdouble)
    ccall((:gsw_cp_ice, libgswteos), Cdouble, (Cdouble, Cdouble), t, p)
end

function gsw_cp_t_exact(sa::Cdouble, t::Cdouble, p::Cdouble)
    ccall((:gsw_cp_t_exact, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, t, p)
end

function gsw_ct_first_derivatives(sa::Cdouble, pt::Cdouble, ct_sa, ct_pt)
    ccall((:gsw_ct_first_derivatives, libgswteos), Nothing, (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), sa, pt, ct_sa, ct_pt)
end

function gsw_ct_first_derivatives_wrt_t_exact(sa::Cdouble, t::Cdouble, p::Cdouble, ct_sa_wrt_t, ct_t_wrt_t, ct_p_wrt_t)
    ccall((:gsw_ct_first_derivatives_wrt_t_exact, libgswteos), Nothing, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, t, p, ct_sa_wrt_t, ct_t_wrt_t, ct_p_wrt_t)
end

function gsw_ct_freezing(sa::Cdouble, p::Cdouble, saturation_fraction::Cdouble)
    ccall((:gsw_ct_freezing, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, p, saturation_fraction)
end

function gsw_ct_freezing_first_derivatives(sa::Cdouble, p::Cdouble, saturation_fraction::Cdouble, ctfreezing_sa, ctfreezing_p)
    ccall((:gsw_ct_freezing_first_derivatives, libgswteos), Nothing, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), sa, p, saturation_fraction, ctfreezing_sa, ctfreezing_p)
end

function gsw_ct_freezing_first_derivatives_poly(sa::Cdouble, p::Cdouble, saturation_fraction::Cdouble, ctfreezing_sa, ctfreezing_p)
    ccall((:gsw_ct_freezing_first_derivatives_poly, libgswteos), Nothing, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), sa, p, saturation_fraction, ctfreezing_sa, ctfreezing_p)
end

function gsw_ct_freezing_poly(sa::Cdouble, p::Cdouble, saturation_fraction::Cdouble)
    ccall((:gsw_ct_freezing_poly, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, p, saturation_fraction)
end

function gsw_ct_from_enthalpy(sa::Cdouble, h::Cdouble, p::Cdouble)
    ccall((:gsw_ct_from_enthalpy, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, h, p)
end

function gsw_ct_from_enthalpy_exact(sa::Cdouble, h::Cdouble, p::Cdouble)
    ccall((:gsw_ct_from_enthalpy_exact, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, h, p)
end

function gsw_ct_from_entropy(sa::Cdouble, entropy::Cdouble)
    ccall((:gsw_ct_from_entropy, libgswteos), Cdouble, (Cdouble, Cdouble), sa, entropy)
end

function gsw_ct_from_pt(sa::Cdouble, pt::Cdouble)
    ccall((:gsw_ct_from_pt, libgswteos), Cdouble, (Cdouble, Cdouble), sa, pt)
end

function gsw_ct_from_rho(rho::Cdouble, sa::Cdouble, p::Cdouble, ct, ct_multiple)
    ccall((:gsw_ct_from_rho, libgswteos), Nothing, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), rho, sa, p, ct, ct_multiple)
end

function gsw_ct_from_t(sa::Cdouble, t::Cdouble, p::Cdouble)
    ccall((:gsw_ct_from_t, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, t, p)
end

function gsw_ct_maxdensity(sa::Cdouble, p::Cdouble)
    ccall((:gsw_ct_maxdensity, libgswteos), Cdouble, (Cdouble, Cdouble), sa, p)
end

function gsw_ct_second_derivatives(sa::Cdouble, pt::Cdouble, ct_sa_sa, ct_sa_pt, ct_pt_pt)
    ccall((:gsw_ct_second_derivatives, libgswteos), Nothing, (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, pt, ct_sa_sa, ct_sa_pt, ct_pt_pt)
end

function gsw_deltasa_atlas(p::Cdouble, lon::Cdouble, lat::Cdouble)
    ccall((:gsw_deltasa_atlas, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), p, lon, lat)
end

function gsw_deltasa_from_sp(sp::Cdouble, p::Cdouble, lon::Cdouble, lat::Cdouble)
    ccall((:gsw_deltasa_from_sp, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble, Cdouble), sp, p, lon, lat)
end

function gsw_dilution_coefficient_t_exact(sa::Cdouble, t::Cdouble, p::Cdouble)
    ccall((:gsw_dilution_coefficient_t_exact, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, t, p)
end

function gsw_dynamic_enthalpy(sa::Cdouble, ct::Cdouble, p::Cdouble)
    ccall((:gsw_dynamic_enthalpy, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, ct, p)
end

function gsw_enthalpy_ct_exact(sa::Cdouble, ct::Cdouble, p::Cdouble)
    ccall((:gsw_enthalpy_ct_exact, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, ct, p)
end

function gsw_enthalpy_diff(sa::Cdouble, ct::Cdouble, p_shallow::Cdouble, p_deep::Cdouble)
    ccall((:gsw_enthalpy_diff, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble, Cdouble), sa, ct, p_shallow, p_deep)
end

function gsw_enthalpy(sa::Cdouble, ct::Cdouble, p::Cdouble)
    ccall((:gsw_enthalpy, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, ct, p)
end

function gsw_enthalpy_first_derivatives_ct_exact(sa::Cdouble, ct::Cdouble, p::Cdouble, h_sa, h_ct)
    ccall((:gsw_enthalpy_first_derivatives_ct_exact, libgswteos), Nothing, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, h_sa, h_ct)
end

function gsw_enthalpy_first_derivatives(sa::Cdouble, ct::Cdouble, p::Cdouble, h_sa, h_ct)
    ccall((:gsw_enthalpy_first_derivatives, libgswteos), Nothing, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, h_sa, h_ct)
end

function gsw_enthalpy_ice(t::Cdouble, p::Cdouble)
    ccall((:gsw_enthalpy_ice, libgswteos), Cdouble, (Cdouble, Cdouble), t, p)
end

function gsw_enthalpy_second_derivatives_ct_exact(sa::Cdouble, ct::Cdouble, p::Cdouble, h_sa_sa, h_sa_ct, h_ct_ct)
    ccall((:gsw_enthalpy_second_derivatives_ct_exact, libgswteos), Nothing, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, h_sa_sa, h_sa_ct, h_ct_ct)
end

function gsw_enthalpy_second_derivatives(sa::Cdouble, ct::Cdouble, p::Cdouble, h_sa_sa, h_sa_ct, h_ct_ct)
    ccall((:gsw_enthalpy_second_derivatives, libgswteos), Nothing, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, h_sa_sa, h_sa_ct, h_ct_ct)
end

function gsw_enthalpy_sso_0(p::Cdouble)
    ccall((:gsw_enthalpy_sso_0, libgswteos), Cdouble, (Cdouble,), p)
end

function gsw_enthalpy_t_exact(sa::Cdouble, t::Cdouble, p::Cdouble)
    ccall((:gsw_enthalpy_t_exact, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, t, p)
end

function gsw_entropy_first_derivatives(sa::Cdouble, ct::Cdouble, eta_sa, eta_ct)
    ccall((:gsw_entropy_first_derivatives, libgswteos), Nothing, (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, eta_sa, eta_ct)
end

function gsw_entropy_from_ct(sa::Cdouble, ct::Cdouble)
    ccall((:gsw_entropy_from_ct, libgswteos), Cdouble, (Cdouble, Cdouble), sa, ct)
end

function gsw_entropy_from_pt(sa::Cdouble, pt::Cdouble)
    ccall((:gsw_entropy_from_pt, libgswteos), Cdouble, (Cdouble, Cdouble), sa, pt)
end

function gsw_entropy_from_t(sa::Cdouble, t::Cdouble, p::Cdouble)
    ccall((:gsw_entropy_from_t, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, t, p)
end

function gsw_entropy_ice(t::Cdouble, p::Cdouble)
    ccall((:gsw_entropy_ice, libgswteos), Cdouble, (Cdouble, Cdouble), t, p)
end

function gsw_entropy_part(sa::Cdouble, t::Cdouble, p::Cdouble)
    ccall((:gsw_entropy_part, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, t, p)
end

function gsw_entropy_part_zerop(sa::Cdouble, pt0::Cdouble)
    ccall((:gsw_entropy_part_zerop, libgswteos), Cdouble, (Cdouble, Cdouble), sa, pt0)
end

function gsw_entropy_second_derivatives(sa::Cdouble, ct::Cdouble, eta_sa_sa, eta_sa_ct, eta_ct_ct)
    ccall((:gsw_entropy_second_derivatives, libgswteos), Nothing, (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, eta_sa_sa, eta_sa_ct, eta_ct_ct)
end

function gsw_fdelta(p::Cdouble, lon::Cdouble, lat::Cdouble)
    ccall((:gsw_fdelta, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), p, lon, lat)
end

function gsw_frazil_properties(sa_bulk::Cdouble, h_bulk::Cdouble, p::Cdouble, sa_final, ct_final, w_ih_final)
    ccall((:gsw_frazil_properties, libgswteos), Nothing, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa_bulk, h_bulk, p, sa_final, ct_final, w_ih_final)
end

function gsw_frazil_properties_potential(sa_bulk::Cdouble, h_pot_bulk::Cdouble, p::Cdouble, sa_final, ct_final, w_ih_final)
    ccall((:gsw_frazil_properties_potential, libgswteos), Nothing, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa_bulk, h_pot_bulk, p, sa_final, ct_final, w_ih_final)
end

function gsw_frazil_properties_potential_poly(sa_bulk::Cdouble, h_pot_bulk::Cdouble, p::Cdouble, sa_final, ct_final, w_ih_final)
    ccall((:gsw_frazil_properties_potential_poly, libgswteos), Nothing, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa_bulk, h_pot_bulk, p, sa_final, ct_final, w_ih_final)
end

function gsw_frazil_ratios_adiabatic(sa::Cdouble, p::Cdouble, w_ih::Cdouble, dsa_dct_frazil, dsa_dp_frazil, dct_dp_frazil)
    ccall((:gsw_frazil_ratios_adiabatic, libgswteos), Nothing, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, p, w_ih, dsa_dct_frazil, dsa_dp_frazil, dct_dp_frazil)
end

function gsw_frazil_ratios_adiabatic_poly(sa::Cdouble, p::Cdouble, w_ih::Cdouble, dsa_dct_frazil, dsa_dp_frazil, dct_dp_frazil)
    ccall((:gsw_frazil_ratios_adiabatic_poly, libgswteos), Nothing, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, p, w_ih, dsa_dct_frazil, dsa_dp_frazil, dct_dp_frazil)
end

function gsw_geo_strf_dyn_height(sa, ct, p, p_ref::Cdouble, n_levels::Cint, dyn_height)
    ccall((:gsw_geo_strf_dyn_height, libgswteos), Ptr{Cdouble}, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cint, Ptr{Cdouble}), sa, ct, p, p_ref, n_levels, dyn_height)
end

function gsw_geo_strf_dyn_height_1(sa, ct, p, p_ref::Cdouble, n_levels::Cint, dyn_height, max_dp_i::Cdouble, interp_method::Cint)
    ccall((:gsw_geo_strf_dyn_height_1, libgswteos), Cint, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cint, Ptr{Cdouble}, Cdouble, Cint), sa, ct, p, p_ref, n_levels, dyn_height, max_dp_i, interp_method)
end

function gsw_geo_strf_dyn_height_pc(sa, ct, delta_p, n_levels::Cint, geo_strf_dyn_height_pc, p_mid)
    ccall((:gsw_geo_strf_dyn_height_pc, libgswteos), Ptr{Cdouble}, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, delta_p, n_levels, geo_strf_dyn_height_pc, p_mid)
end

function gsw_gibbs_ice(nt::Cint, np::Cint, t::Cdouble, p::Cdouble)
    ccall((:gsw_gibbs_ice, libgswteos), Cdouble, (Cint, Cint, Cdouble, Cdouble), nt, np, t, p)
end

function gsw_gibbs_ice_part_t(t::Cdouble, p::Cdouble)
    ccall((:gsw_gibbs_ice_part_t, libgswteos), Cdouble, (Cdouble, Cdouble), t, p)
end

function gsw_gibbs_ice_pt0(pt0::Cdouble)
    ccall((:gsw_gibbs_ice_pt0, libgswteos), Cdouble, (Cdouble,), pt0)
end

function gsw_gibbs_ice_pt0_pt0(pt0::Cdouble)
    ccall((:gsw_gibbs_ice_pt0_pt0, libgswteos), Cdouble, (Cdouble,), pt0)
end

function gsw_gibbs(ns::Cint, nt::Cint, np::Cint, sa::Cdouble, t::Cdouble, p::Cdouble)
    ccall((:gsw_gibbs, libgswteos), Cdouble, (Cint, Cint, Cint, Cdouble, Cdouble, Cdouble), ns, nt, np, sa, t, p)
end

function gsw_gibbs_pt0_pt0(sa::Cdouble, pt0::Cdouble)
    ccall((:gsw_gibbs_pt0_pt0, libgswteos), Cdouble, (Cdouble, Cdouble), sa, pt0)
end

function gsw_grav(lat::Cdouble, p::Cdouble)
    ccall((:gsw_grav, libgswteos), Cdouble, (Cdouble, Cdouble), lat, p)
end

function gsw_helmholtz_energy_ice(t::Cdouble, p::Cdouble)
    ccall((:gsw_helmholtz_energy_ice, libgswteos), Cdouble, (Cdouble, Cdouble), t, p)
end

function gsw_hill_ratio_at_sp2(t::Cdouble)
    ccall((:gsw_hill_ratio_at_sp2, libgswteos), Cdouble, (Cdouble,), t)
end

function gsw_ice_fraction_to_freeze_seawater(sa::Cdouble, ct::Cdouble, p::Cdouble, t_ih::Cdouble, sa_freeze, ct_freeze, w_ih)
    ccall((:gsw_ice_fraction_to_freeze_seawater, libgswteos), Nothing, (Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, t_ih, sa_freeze, ct_freeze, w_ih)
end

function gsw_internal_energy(sa::Cdouble, ct::Cdouble, p::Cdouble)
    ccall((:gsw_internal_energy, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, ct, p)
end

function gsw_internal_energy_ice(t::Cdouble, p::Cdouble)
    ccall((:gsw_internal_energy_ice, libgswteos), Cdouble, (Cdouble, Cdouble), t, p)
end

function gsw_ipv_vs_fnsquared_ratio(sa, ct, p, p_ref::Cdouble, nz::Cint, ipv_vs_fnsquared_ratio, p_mid)
    ccall((:gsw_ipv_vs_fnsquared_ratio, libgswteos), Nothing, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cint, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, p_ref, nz, ipv_vs_fnsquared_ratio, p_mid)
end

function gsw_kappa_const_t_ice(t::Cdouble, p::Cdouble)
    ccall((:gsw_kappa_const_t_ice, libgswteos), Cdouble, (Cdouble, Cdouble), t, p)
end

function gsw_kappa(sa::Cdouble, ct::Cdouble, p::Cdouble)
    ccall((:gsw_kappa, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, ct, p)
end

function gsw_kappa_ice(t::Cdouble, p::Cdouble)
    ccall((:gsw_kappa_ice, libgswteos), Cdouble, (Cdouble, Cdouble), t, p)
end

function gsw_kappa_t_exact(sa::Cdouble, t::Cdouble, p::Cdouble)
    ccall((:gsw_kappa_t_exact, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, t, p)
end

function gsw_latentheat_evap_ct(sa::Cdouble, ct::Cdouble)
    ccall((:gsw_latentheat_evap_ct, libgswteos), Cdouble, (Cdouble, Cdouble), sa, ct)
end

function gsw_latentheat_evap_t(sa::Cdouble, t::Cdouble)
    ccall((:gsw_latentheat_evap_t, libgswteos), Cdouble, (Cdouble, Cdouble), sa, t)
end

function gsw_latentheat_melting(sa::Cdouble, p::Cdouble)
    ccall((:gsw_latentheat_melting, libgswteos), Cdouble, (Cdouble, Cdouble), sa, p)
end

function gsw_linear_interp_sa_ct(sa, ct, p, np::Cint, p_i, npi::Cint, sa_i, ct_i)
    ccall((:gsw_linear_interp_sa_ct, libgswteos), Nothing, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, np, p_i, npi, sa_i, ct_i)
end

function gsw_melting_ice_equilibrium_sa_ct_ratio(sa::Cdouble, p::Cdouble)
    ccall((:gsw_melting_ice_equilibrium_sa_ct_ratio, libgswteos), Cdouble, (Cdouble, Cdouble), sa, p)
end

function gsw_melting_ice_equilibrium_sa_ct_ratio_poly(sa::Cdouble, p::Cdouble)
    ccall((:gsw_melting_ice_equilibrium_sa_ct_ratio_poly, libgswteos), Cdouble, (Cdouble, Cdouble), sa, p)
end

function gsw_melting_ice_into_seawater(sa::Cdouble, ct::Cdouble, p::Cdouble, w_ih::Cdouble, t_ih::Cdouble, sa_final, ct_final, w_ih_final)
    ccall((:gsw_melting_ice_into_seawater, libgswteos), Nothing, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, w_ih, t_ih, sa_final, ct_final, w_ih_final)
end

function gsw_melting_ice_sa_ct_ratio(sa::Cdouble, ct::Cdouble, p::Cdouble, t_ih::Cdouble)
    ccall((:gsw_melting_ice_sa_ct_ratio, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble, Cdouble), sa, ct, p, t_ih)
end

function gsw_melting_ice_sa_ct_ratio_poly(sa::Cdouble, ct::Cdouble, p::Cdouble, t_ih::Cdouble)
    ccall((:gsw_melting_ice_sa_ct_ratio_poly, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble, Cdouble), sa, ct, p, t_ih)
end

function gsw_melting_seaice_equilibrium_sa_ct_ratio(sa::Cdouble, p::Cdouble)
    ccall((:gsw_melting_seaice_equilibrium_sa_ct_ratio, libgswteos), Cdouble, (Cdouble, Cdouble), sa, p)
end

function gsw_melting_seaice_equilibrium_sa_ct_ratio_poly(sa::Cdouble, p::Cdouble)
    ccall((:gsw_melting_seaice_equilibrium_sa_ct_ratio_poly, libgswteos), Cdouble, (Cdouble, Cdouble), sa, p)
end

function gsw_melting_seaice_into_seawater(sa::Cdouble, ct::Cdouble, p::Cdouble, w_seaice::Cdouble, sa_seaice::Cdouble, t_seaice::Cdouble, sa_final, ct_final)
    ccall((:gsw_melting_seaice_into_seawater, libgswteos), Nothing, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, w_seaice, sa_seaice, t_seaice, sa_final, ct_final)
end

function gsw_melting_seaice_sa_ct_ratio(sa::Cdouble, ct::Cdouble, p::Cdouble, sa_seaice::Cdouble, t_seaice::Cdouble)
    ccall((:gsw_melting_seaice_sa_ct_ratio, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble), sa, ct, p, sa_seaice, t_seaice)
end

function gsw_melting_seaice_sa_ct_ratio_poly(sa::Cdouble, ct::Cdouble, p::Cdouble, sa_seaice::Cdouble, t_seaice::Cdouble)
    ccall((:gsw_melting_seaice_sa_ct_ratio_poly, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble), sa, ct, p, sa_seaice, t_seaice)
end

function gsw_nsquared(sa, ct, p, lat, nz::Cint, n2, p_mid)
    ccall((:gsw_nsquared, libgswteos), Nothing, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, lat, nz, n2, p_mid)
end

function gsw_o2sol(sa::Cdouble, ct::Cdouble, p::Cdouble, lon::Cdouble, lat::Cdouble)
    ccall((:gsw_o2sol, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble), sa, ct, p, lon, lat)
end

function gsw_o2sol_sp_pt(sp::Cdouble, pt::Cdouble)
    ccall((:gsw_o2sol_sp_pt, libgswteos), Cdouble, (Cdouble, Cdouble), sp, pt)
end

function gsw_pot_enthalpy_from_pt_ice(pt0_ice::Cdouble)
    ccall((:gsw_pot_enthalpy_from_pt_ice, libgswteos), Cdouble, (Cdouble,), pt0_ice)
end

function gsw_pot_enthalpy_from_pt_ice_poly(pt0_ice::Cdouble)
    ccall((:gsw_pot_enthalpy_from_pt_ice_poly, libgswteos), Cdouble, (Cdouble,), pt0_ice)
end

function gsw_pot_enthalpy_ice_freezing(sa::Cdouble, p::Cdouble)
    ccall((:gsw_pot_enthalpy_ice_freezing, libgswteos), Cdouble, (Cdouble, Cdouble), sa, p)
end

function gsw_pot_enthalpy_ice_freezing_first_derivatives(sa::Cdouble, p::Cdouble, pot_enthalpy_ice_freezing_sa, pot_enthalpy_ice_freezing_p)
    ccall((:gsw_pot_enthalpy_ice_freezing_first_derivatives, libgswteos), Nothing, (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), sa, p, pot_enthalpy_ice_freezing_sa, pot_enthalpy_ice_freezing_p)
end

function gsw_pot_enthalpy_ice_freezing_first_derivatives_poly(sa::Cdouble, p::Cdouble, pot_enthalpy_ice_freezing_sa, pot_enthalpy_ice_freezing_p)
    ccall((:gsw_pot_enthalpy_ice_freezing_first_derivatives_poly, libgswteos), Nothing, (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), sa, p, pot_enthalpy_ice_freezing_sa, pot_enthalpy_ice_freezing_p)
end

function gsw_pot_enthalpy_ice_freezing_poly(sa::Cdouble, p::Cdouble)
    ccall((:gsw_pot_enthalpy_ice_freezing_poly, libgswteos), Cdouble, (Cdouble, Cdouble), sa, p)
end

function gsw_pot_rho_t_exact(sa::Cdouble, t::Cdouble, p::Cdouble, p_ref::Cdouble)
    ccall((:gsw_pot_rho_t_exact, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble, Cdouble), sa, t, p, p_ref)
end

function gsw_pressure_coefficient_ice(t::Cdouble, p::Cdouble)
    ccall((:gsw_pressure_coefficient_ice, libgswteos), Cdouble, (Cdouble, Cdouble), t, p)
end

function gsw_pressure_freezing_ct(sa::Cdouble, ct::Cdouble, saturation_fraction::Cdouble)
    ccall((:gsw_pressure_freezing_ct, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, ct, saturation_fraction)
end

function gsw_pt0_cold_ice_poly(pot_enthalpy_ice::Cdouble)
    ccall((:gsw_pt0_cold_ice_poly, libgswteos), Cdouble, (Cdouble,), pot_enthalpy_ice)
end

function gsw_pt0_from_t(sa::Cdouble, t::Cdouble, p::Cdouble)
    ccall((:gsw_pt0_from_t, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, t, p)
end

function gsw_pt0_from_t_ice(t::Cdouble, p::Cdouble)
    ccall((:gsw_pt0_from_t_ice, libgswteos), Cdouble, (Cdouble, Cdouble), t, p)
end

function gsw_pt_first_derivatives(sa::Cdouble, ct::Cdouble, pt_sa, pt_ct)
    ccall((:gsw_pt_first_derivatives, libgswteos), Nothing, (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, pt_sa, pt_ct)
end

function gsw_pt_from_ct(sa::Cdouble, ct::Cdouble)
    ccall((:gsw_pt_from_ct, libgswteos), Cdouble, (Cdouble, Cdouble), sa, ct)
end

function gsw_pt_from_entropy(sa::Cdouble, entropy::Cdouble)
    ccall((:gsw_pt_from_entropy, libgswteos), Cdouble, (Cdouble, Cdouble), sa, entropy)
end

function gsw_pt_from_pot_enthalpy_ice(pot_enthalpy_ice::Cdouble)
    ccall((:gsw_pt_from_pot_enthalpy_ice, libgswteos), Cdouble, (Cdouble,), pot_enthalpy_ice)
end

function gsw_pt_from_pot_enthalpy_ice_poly_dh(pot_enthalpy_ice::Cdouble)
    ccall((:gsw_pt_from_pot_enthalpy_ice_poly_dh, libgswteos), Cdouble, (Cdouble,), pot_enthalpy_ice)
end

function gsw_pt_from_pot_enthalpy_ice_poly(pot_enthalpy_ice::Cdouble)
    ccall((:gsw_pt_from_pot_enthalpy_ice_poly, libgswteos), Cdouble, (Cdouble,), pot_enthalpy_ice)
end

function gsw_pt_from_t(sa::Cdouble, t::Cdouble, p::Cdouble, p_ref::Cdouble)
    ccall((:gsw_pt_from_t, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble, Cdouble), sa, t, p, p_ref)
end

function gsw_pt_from_t_ice(t::Cdouble, p::Cdouble, p_ref::Cdouble)
    ccall((:gsw_pt_from_t_ice, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), t, p, p_ref)
end

function gsw_pt_second_derivatives(sa::Cdouble, ct::Cdouble, pt_sa_sa, pt_sa_ct, pt_ct_ct)
    ccall((:gsw_pt_second_derivatives, libgswteos), Nothing, (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, pt_sa_sa, pt_sa_ct, pt_ct_ct)
end

function gsw_rho_alpha_beta(sa::Cdouble, ct::Cdouble, p::Cdouble, rho, alpha, beta)
    ccall((:gsw_rho_alpha_beta, libgswteos), Nothing, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, rho, alpha, beta)
end

function gsw_rho(sa::Cdouble, ct::Cdouble, p::Cdouble)
    ccall((:gsw_rho, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, ct, p)
end

function gsw_rho_first_derivatives(sa::Cdouble, ct::Cdouble, p::Cdouble, drho_dsa, drho_dct, drho_dp)
    ccall((:gsw_rho_first_derivatives, libgswteos), Nothing, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, drho_dsa, drho_dct, drho_dp)
end

function gsw_rho_first_derivatives_wrt_enthalpy(sa::Cdouble, ct::Cdouble, p::Cdouble, rho_sa, rho_h)
    ccall((:gsw_rho_first_derivatives_wrt_enthalpy, libgswteos), Nothing, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, rho_sa, rho_h)
end

function gsw_rho_ice(t::Cdouble, p::Cdouble)
    ccall((:gsw_rho_ice, libgswteos), Cdouble, (Cdouble, Cdouble), t, p)
end

function gsw_rho_second_derivatives(sa::Cdouble, ct::Cdouble, p::Cdouble, rho_sa_sa, rho_sa_ct, rho_ct_ct, rho_sa_p, rho_ct_p)
    ccall((:gsw_rho_second_derivatives, libgswteos), Nothing, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, rho_sa_sa, rho_sa_ct, rho_ct_ct, rho_sa_p, rho_ct_p)
end

function gsw_rho_second_derivatives_wrt_enthalpy(sa::Cdouble, ct::Cdouble, p::Cdouble, rho_sa_sa, rho_sa_h, rho_h_h)
    ccall((:gsw_rho_second_derivatives_wrt_enthalpy, libgswteos), Nothing, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, rho_sa_sa, rho_sa_h, rho_h_h)
end

function gsw_rho_t_exact(sa::Cdouble, t::Cdouble, p::Cdouble)
    ccall((:gsw_rho_t_exact, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, t, p)
end

function gsw_rr68_interp_sa_ct(sa, ct, p, mp::Cint, p_i, mp_i::Cint, sa_i, ct_i)
    ccall((:gsw_rr68_interp_sa_ct, libgswteos), Nothing, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, mp, p_i, mp_i, sa_i, ct_i)
end

function gsw_saar(p::Cdouble, lon::Cdouble, lat::Cdouble)
    ccall((:gsw_saar, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), p, lon, lat)
end

function gsw_sa_freezing_estimate(p::Cdouble, saturation_fraction::Cdouble, ct, t)
    ccall((:gsw_sa_freezing_estimate, libgswteos), Cdouble, (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), p, saturation_fraction, ct, t)
end

function gsw_sa_freezing_from_ct(ct::Cdouble, p::Cdouble, saturation_fraction::Cdouble)
    ccall((:gsw_sa_freezing_from_ct, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), ct, p, saturation_fraction)
end

function gsw_sa_freezing_from_ct_poly(ct::Cdouble, p::Cdouble, saturation_fraction::Cdouble)
    ccall((:gsw_sa_freezing_from_ct_poly, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), ct, p, saturation_fraction)
end

function gsw_sa_freezing_from_t(t::Cdouble, p::Cdouble, saturation_fraction::Cdouble)
    ccall((:gsw_sa_freezing_from_t, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), t, p, saturation_fraction)
end

function gsw_sa_freezing_from_t_poly(t::Cdouble, p::Cdouble, saturation_fraction::Cdouble)
    ccall((:gsw_sa_freezing_from_t_poly, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), t, p, saturation_fraction)
end

function gsw_sa_from_rho(rho::Cdouble, ct::Cdouble, p::Cdouble)
    ccall((:gsw_sa_from_rho, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), rho, ct, p)
end

function gsw_sa_from_sp_baltic(sp::Cdouble, lon::Cdouble, lat::Cdouble)
    ccall((:gsw_sa_from_sp_baltic, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sp, lon, lat)
end

function gsw_sa_from_sp(sp::Cdouble, p::Cdouble, lon::Cdouble, lat::Cdouble)
    ccall((:gsw_sa_from_sp, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble, Cdouble), sp, p, lon, lat)
end

function gsw_sa_from_sstar(sstar::Cdouble, p::Cdouble, lon::Cdouble, lat::Cdouble)
    ccall((:gsw_sa_from_sstar, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble, Cdouble), sstar, p, lon, lat)
end

function gsw_sa_p_inrange(sa::Cdouble, p::Cdouble)
    ccall((:gsw_sa_p_inrange, libgswteos), Cint, (Cdouble, Cdouble), sa, p)
end

function gsw_seaice_fraction_to_freeze_seawater(sa::Cdouble, ct::Cdouble, p::Cdouble, sa_seaice::Cdouble, t_seaice::Cdouble, sa_freeze, ct_freeze, w_seaice)
    ccall((:gsw_seaice_fraction_to_freeze_seawater, libgswteos), Nothing, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, sa_seaice, t_seaice, sa_freeze, ct_freeze, w_seaice)
end

function gsw_sigma0(sa::Cdouble, ct::Cdouble)
    ccall((:gsw_sigma0, libgswteos), Cdouble, (Cdouble, Cdouble), sa, ct)
end

function gsw_sigma1(sa::Cdouble, ct::Cdouble)
    ccall((:gsw_sigma1, libgswteos), Cdouble, (Cdouble, Cdouble), sa, ct)
end

function gsw_sigma2(sa::Cdouble, ct::Cdouble)
    ccall((:gsw_sigma2, libgswteos), Cdouble, (Cdouble, Cdouble), sa, ct)
end

function gsw_sigma3(sa::Cdouble, ct::Cdouble)
    ccall((:gsw_sigma3, libgswteos), Cdouble, (Cdouble, Cdouble), sa, ct)
end

function gsw_sigma4(sa::Cdouble, ct::Cdouble)
    ccall((:gsw_sigma4, libgswteos), Cdouble, (Cdouble, Cdouble), sa, ct)
end

function gsw_sound_speed(sa::Cdouble, ct::Cdouble, p::Cdouble)
    ccall((:gsw_sound_speed, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, ct, p)
end

function gsw_sound_speed_ice(t::Cdouble, p::Cdouble)
    ccall((:gsw_sound_speed_ice, libgswteos), Cdouble, (Cdouble, Cdouble), t, p)
end

function gsw_sound_speed_t_exact(sa::Cdouble, t::Cdouble, p::Cdouble)
    ccall((:gsw_sound_speed_t_exact, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, t, p)
end

function gsw_specvol_alpha_beta(sa::Cdouble, ct::Cdouble, p::Cdouble, specvol, alpha, beta)
    ccall((:gsw_specvol_alpha_beta, libgswteos), Nothing, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, specvol, alpha, beta)
end

function gsw_specvol_anom_standard(sa::Cdouble, ct::Cdouble, p::Cdouble)
    ccall((:gsw_specvol_anom_standard, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, ct, p)
end

function gsw_specvol(sa::Cdouble, ct::Cdouble, p::Cdouble)
    ccall((:gsw_specvol, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, ct, p)
end

function gsw_specvol_first_derivatives(sa::Cdouble, ct::Cdouble, p::Cdouble, v_sa, v_ct, v_p)
    ccall((:gsw_specvol_first_derivatives, libgswteos), Nothing, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, v_sa, v_ct, v_p)
end

function gsw_specvol_first_derivatives_wrt_enthalpy(sa::Cdouble, ct::Cdouble, p::Cdouble, v_sa, v_h)
    ccall((:gsw_specvol_first_derivatives_wrt_enthalpy, libgswteos), Nothing, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, v_sa, v_h)
end

function gsw_specvol_ice(t::Cdouble, p::Cdouble)
    ccall((:gsw_specvol_ice, libgswteos), Cdouble, (Cdouble, Cdouble), t, p)
end

function gsw_specvol_second_derivatives(sa::Cdouble, ct::Cdouble, p::Cdouble, v_sa_sa, v_sa_ct, v_ct_ct, v_sa_p, v_ct_p)
    ccall((:gsw_specvol_second_derivatives, libgswteos), Nothing, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, v_sa_sa, v_sa_ct, v_ct_ct, v_sa_p, v_ct_p)
end

function gsw_specvol_second_derivatives_wrt_enthalpy(sa::Cdouble, ct::Cdouble, p::Cdouble, v_sa_sa, v_sa_h, v_h_h)
    ccall((:gsw_specvol_second_derivatives_wrt_enthalpy, libgswteos), Nothing, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, v_sa_sa, v_sa_h, v_h_h)
end

function gsw_specvol_sso_0(p::Cdouble)
    ccall((:gsw_specvol_sso_0, libgswteos), Cdouble, (Cdouble,), p)
end

function gsw_specvol_t_exact(sa::Cdouble, t::Cdouble, p::Cdouble)
    ccall((:gsw_specvol_t_exact, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, t, p)
end

function gsw_sp_from_c(c::Cdouble, t::Cdouble, p::Cdouble)
    ccall((:gsw_sp_from_c, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), c, t, p)
end

function gsw_sp_from_sa_baltic(sa::Cdouble, lon::Cdouble, lat::Cdouble)
    ccall((:gsw_sp_from_sa_baltic, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, lon, lat)
end

function gsw_sp_from_sa(sa::Cdouble, p::Cdouble, lon::Cdouble, lat::Cdouble)
    ccall((:gsw_sp_from_sa, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble, Cdouble), sa, p, lon, lat)
end

function gsw_sp_from_sk(sk::Cdouble)
    ccall((:gsw_sp_from_sk, libgswteos), Cdouble, (Cdouble,), sk)
end

function gsw_sp_from_sr(sr::Cdouble)
    ccall((:gsw_sp_from_sr, libgswteos), Cdouble, (Cdouble,), sr)
end

function gsw_sp_from_sstar(sstar::Cdouble, p::Cdouble, lon::Cdouble, lat::Cdouble)
    ccall((:gsw_sp_from_sstar, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble, Cdouble), sstar, p, lon, lat)
end

function gsw_spiciness0(sa::Cdouble, ct::Cdouble)
    ccall((:gsw_spiciness0, libgswteos), Cdouble, (Cdouble, Cdouble), sa, ct)
end

function gsw_spiciness1(sa::Cdouble, ct::Cdouble)
    ccall((:gsw_spiciness1, libgswteos), Cdouble, (Cdouble, Cdouble), sa, ct)
end

function gsw_spiciness2(sa::Cdouble, ct::Cdouble)
    ccall((:gsw_spiciness2, libgswteos), Cdouble, (Cdouble, Cdouble), sa, ct)
end

function gsw_sr_from_sp(sp::Cdouble)
    ccall((:gsw_sr_from_sp, libgswteos), Cdouble, (Cdouble,), sp)
end

function gsw_sstar_from_sa(sa::Cdouble, p::Cdouble, lon::Cdouble, lat::Cdouble)
    ccall((:gsw_sstar_from_sa, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble, Cdouble), sa, p, lon, lat)
end

function gsw_sstar_from_sp(sp::Cdouble, p::Cdouble, lon::Cdouble, lat::Cdouble)
    ccall((:gsw_sstar_from_sp, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble, Cdouble), sp, p, lon, lat)
end

function gsw_t_deriv_chem_potential_water_t_exact(sa::Cdouble, t::Cdouble, p::Cdouble)
    ccall((:gsw_t_deriv_chem_potential_water_t_exact, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, t, p)
end

function gsw_t_freezing(sa::Cdouble, p::Cdouble, saturation_fraction::Cdouble)
    ccall((:gsw_t_freezing, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, p, saturation_fraction)
end

function gsw_t_freezing_first_derivatives_poly(sa::Cdouble, p::Cdouble, saturation_fraction::Cdouble, tfreezing_sa, tfreezing_p)
    ccall((:gsw_t_freezing_first_derivatives_poly, libgswteos), Nothing, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), sa, p, saturation_fraction, tfreezing_sa, tfreezing_p)
end

function gsw_t_freezing_first_derivatives(sa::Cdouble, p::Cdouble, saturation_fraction::Cdouble, tfreezing_sa, tfreezing_p)
    ccall((:gsw_t_freezing_first_derivatives, libgswteos), Nothing, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), sa, p, saturation_fraction, tfreezing_sa, tfreezing_p)
end

function gsw_t_freezing_poly(sa::Cdouble, p::Cdouble, saturation_fraction::Cdouble)
    ccall((:gsw_t_freezing_poly, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, p, saturation_fraction)
end

function gsw_t_from_ct(sa::Cdouble, ct::Cdouble, p::Cdouble)
    ccall((:gsw_t_from_ct, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, ct, p)
end

function gsw_t_from_pt0_ice(pt0_ice::Cdouble, p::Cdouble)
    ccall((:gsw_t_from_pt0_ice, libgswteos), Cdouble, (Cdouble, Cdouble), pt0_ice, p)
end

function gsw_thermobaric(sa::Cdouble, ct::Cdouble, p::Cdouble)
    ccall((:gsw_thermobaric, libgswteos), Cdouble, (Cdouble, Cdouble, Cdouble), sa, ct, p)
end

function gsw_turner_rsubrho(sa, ct, p, nz::Cint, tu, rsubrho, p_mid)
    ccall((:gsw_turner_rsubrho, libgswteos), Nothing, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, nz, tu, rsubrho, p_mid)
end

function gsw_util_indx(x, n::Cint, z::Cdouble)
    ccall((:gsw_util_indx, libgswteos), Cint, (Ptr{Cdouble}, Cint, Cdouble), x, n, z)
end

function gsw_util_interp1q_int(nx::Cint, x, iy, nxi::Cint, x_i, y_i)
    ccall((:gsw_util_interp1q_int, libgswteos), Ptr{Cdouble}, (Cint, Ptr{Cdouble}, Ptr{Cint}, Cint, Ptr{Cdouble}, Ptr{Cdouble}), nx, x, iy, nxi, x_i, y_i)
end

function gsw_util_linear_interp(nx::Cint, x, ny::Cint, y, nxi::Cint, x_i, y_i)
    ccall((:gsw_util_linear_interp, libgswteos), Ptr{Cdouble}, (Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}), nx, x, ny, y, nxi, x_i, y_i)
end

function gsw_util_sort_real(rarray, nx::Cint, iarray)
    ccall((:gsw_util_sort_real, libgswteos), Nothing, (Ptr{Cdouble}, Cint, Ptr{Cint}), rarray, nx, iarray)
end

function gsw_util_xinterp1(x, y, n::Cint, x0::Cdouble)
    ccall((:gsw_util_xinterp1, libgswteos), Cdouble, (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble), x, y, n, x0)
end

function gsw_util_pchip_interp(x, y, n::Cint, xi, yi, ni::Cint)
    ccall((:gsw_util_pchip_interp, libgswteos), Cint, (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint), x, y, n, xi, yi, ni)
end

function gsw_z_from_p(p::Cdouble, lat::Cdouble)
    ccall((:gsw_z_from_p, libgswteos), Cdouble, (Cdouble, Cdouble), p, lat)
end

function gsw_p_from_z(z::Cdouble, lat::Cdouble)
    ccall((:gsw_p_from_z, libgswteos), Cdouble, (Cdouble, Cdouble), z, lat)
end
