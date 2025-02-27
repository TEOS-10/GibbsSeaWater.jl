function gsw_nsquared(sa::AbstractArray{T}, ct::AbstractArray{T}, p::AbstractArray{T}, lat::AbstractArray{T}) where {T<:Real}
    n_levels = length(sa)
    _res = Array{T}(undef, n_levels - 1)
    _pres = Array{T}(undef, n_levels - 1)
    gsw_nsquared(sa, ct, p, lat, n_levels, _res, _pres)
    return (_res, _pres)
end

function gsw_p_from_z(z, lat)
    return gsw_p_from_z(z, lat, 0.0, 0.0)
end

function gsw_geo_strf_dyn_height(sa::AbstractArray, ct::AbstractArray, p::AbstractArray, p_ref::Real)
    n_levels = length(sa)
    _res = Array{Float64}(undef, n_levels)

    _d = gsw_geo_strf_dyn_height(sa, ct, p, p_ref, n_levels, _res)
    return _res
end

function gsw_geo_strf_dyn_height_1(sa::AbstractArray, ct::AbstractArray, p::AbstractArray, p_ref::Real; max_dp_i=maximum(diff(p)), interp_method=1)
    n_levels = length(sa)
    _res = Array{Float64}(undef, n_levels)

    _d = gsw_geo_strf_dyn_height_1(sa, ct, p, p_ref, n_levels, _res, max_dp_i, interp_method)
    return _res
end

function gsw_geo_strf_dyn_height_pc(sa::AbstractArray, ct::AbstractArray, delta_p::AbstractArray, p_mid::AbstractArray)
    n_levels = length(sa)
    _res = Array{Float64}(undef, n_levels)

    _d = gsw_geo_strf_dyn_height_pc(sa, ct, delta_p, n_levels, _res, p_mid)
    return _res
end

function gsw_turner_rsubrho(sa, ct, p)
    n_levels = length(sa)
    rsubrho = Array{Float64}(undef, n_levels - 1)
    p_mid = Array{Float64}(undef, n_levels - 1)
    tu = Array{Float64}(undef, n_levels - 1)
    ccall((:gsw_turner_rsubrho, libgswteos), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, nz, tu, rsubrho, p_mid)
    return (tu, rsubrho, p_mid)
end

function gsw_add_barrier(input_data, lon, lat, long_grid, lat_grid, dlong_grid, dlat_grid)
    n_levels = 4 #length(input_data)
    output_data = Array{Float64}(undef, n_levels)
    ccall((:gsw_add_barrier, libgswteos), Cvoid, (Ptr{Cdouble}, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}), input_data, lon, lat, long_grid, lat_grid, dlong_grid, dlat_grid, output_data)
    return output_data
end

function gsw_add_mean(data_in)
    n_levels = 4#length(data_in)
    data_out = Array{Float64}(undef, n_levels)
    ccall((:gsw_add_mean, libgswteos), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}), data_in, data_out)
    return data_out
end

function gsw_ct_first_derivatives(sa, pt)
    ct_pt = Ref{Float64}()
    ct_sa = Ref{Float64}()
    ccall((:gsw_ct_first_derivatives, libgswteos), Cvoid, (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), sa, pt, ct_sa, ct_pt)
    ct_pt = ct_pt[]
    ct_sa = ct_sa[]
    return ct_pt, ct_sa
end

function gsw_ct_first_derivatives_wrt_t_exact(sa, t, p)
    ct_sa_wrt_t = Ref{Float64}()
    ct_t_wrt_t = Ref{Float64}()
    ct_p_wrt_t = Ref{Float64}()
    ccall((:gsw_ct_first_derivatives_wrt_t_exact, libgswteos), Cvoid, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, t, p, ct_sa_wrt_t, ct_t_wrt_t, ct_p_wrt_t)
    ct_sa_wrt_t = ct_sa_wrt_t[]
    ct_t_wrt_t = ct_t_wrt_t[]
    ct_p_wrt_t = ct_p_wrt_t[]
    return ct_sa_wrt_t, ct_t_wrt_t, ct_p_wrt_t
end

function gsw_ct_freezing_first_derivatives(sa, p)
    saturation_fraction = Ref{Float64}()
    ctfreezing_p = Ref{Float64}()
    ccall((:gsw_ct_freezing_first_derivatives, libgswteos), Cvoid, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), sa, p, saturation_fraction, ctfreezing_sa, ctfreezing_p)
    saturation_fraction = saturation_fraction[]
    ctfreezing_p = ctfreezing_p[]
    return saturation_fraction, ctfreezing_p
end

function gsw_ct_freezing_first_derivatives_poly(sa, p)
    saturation_fraction = Ref{Float64}()
    ctfreezing_sa = Ref{Float64}()
    ctfreezing_p = Ref{Float64}()
    ccall((:gsw_ct_freezing_first_derivatives_poly, libgswteos), Cvoid, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), sa, p, saturation_fraction, ctfreezing_sa, ctfreezing_p)
    saturation_fraction = saturation_fraction[]
    ctfreezing_sa = ctfreezing_sa[]
    ctfreezing_p = ctfreezing_p[]
    return saturation_fraction, ctfreezing_sa, ctfreezing_p
end

function gsw_ct_from_rho(rho, sa, p)
    ct = Ref{Float64}()
    ct_multiple = Ref{Float64}()
    ccall((:gsw_ct_from_rho, libgswteos), Cvoid, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), rho, sa, p, ct, ct_multiple)
    ct = ct[]
    ct_multiple = ct_multiple[]
    return ct, ct_multiple
end

function gsw_ct_second_derivatives(sa, pt)
    ct_sa_sa = Ref{Float64}()
    ct_sa_pt = Ref{Float64}()
    ct_pt_pt = Ref{Float64}()
    ccall((:gsw_ct_second_derivatives, libgswteos), Cvoid, (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, pt, ct_sa_sa, ct_sa_pt, ct_pt_pt)
    ct_sa_sa = ct_sa_sa[]
    ct_sa_pt = ct_sa_pt[]
    ct_pt_pt = ct_pt_pt[]
    return ct_sa_sa, ct_sa_pt, ct_pt_pt
end

function gsw_enthalpy_first_derivatives(sa, ct, p)
    h_sa = Ref{Float64}()
    h_ct = Ref{Float64}()
    ccall((:gsw_enthalpy_first_derivatives, libgswteos), Cvoid, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, h_sa, h_ct)
    h_sa = h_sa[]
    h_ct = h_ct[]
    return h_sa, h_ct
end

function gsw_enthalpy_second_derivatives_ct_exact(sa, ct, p)
    h_sa_sa = Ref{Float64}()
    h_sa_ct = Ref{Float64}()
    h_ct_ct = Ref{Float64}()
    ccall((:gsw_enthalpy_second_derivatives_ct_exact, libgswteos), Cvoid, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, h_sa_sa, h_sa_ct, h_ct_ct)
    h_sa_sa = h_sa_sa[]
    h_sa_ct = h_sa_ct[]
    h_ct_ct = h_ct_ct[]
    return h_sa_sa, h_sa_ct, h_ct_ct
end

function gsw_enthalpy_second_derivatives(sa, ct, p)
    h_sa_sa = Ref{Float64}()
    h_sa_ct = Ref{Float64}()
    h_ct_ct = Ref{Float64}()
    ccall((:gsw_enthalpy_second_derivatives, libgswteos), Cvoid, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, h_sa_sa, h_sa_ct, h_ct_ct)
    h_sa_sa = h_sa_sa[]
    h_sa_ct = h_sa_ct[]
    h_ct_ct = h_ct_ct[]
    return h_sa_sa, h_sa_ct, h_ct_ct
end

function gsw_enthalpy_first_derivatives_ct_exact(sa, ct, p)
    h_sa = Ref{Float64}()
    h_ct = Ref{Float64}()
    ccall((:gsw_enthalpy_first_derivatives_ct_exact, libgswteos), Cvoid, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, h_sa, h_ct)
    h_sa = h_sa[]
    h_ct = h_ct[]
    return h_sa, h_ct
end

function gsw_entropy_first_derivatives(sa, ct)
    eta_sa = Ref{Float64}()
    eta_ct = Ref{Float64}()
    ccall((:gsw_entropy_first_derivatives, libgswteos), Cvoid, (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, eta_sa, eta_ct)
    eta_sa = eta_sa[]
    eta_ct = eta_ct[]
    return eta_sa, eta_ct
end

function gsw_entropy_second_derivatives(sa, ct)
    eta_sa_sa = Ref{Float64}()
    eta_sa_ct = Ref{Float64}()
    eta_ct_ct = Ref{Float64}()
    ccall((:gsw_entropy_second_derivatives, libgswteos), Cvoid, (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, eta_sa_sa, eta_sa_ct, eta_ct_ct)
    eta_sa_sa = eta_sa_sa[]
    eta_sa_ct = eta_sa_ct[]
    eta_ct_ct = eta_ct_ct[]
    return eta_sa_sa, eta_sa_ct, eta_ct_ct
end

function gsw_frazil_properties(sa_bulk, h_bulk, p)
    sa_final = Ref{Float64}()
    ct_final = Ref{Float64}()
    w_ih_final = Ref{Float64}()
    ccall((:gsw_frazil_properties, libgswteos), Cvoid, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa_bulk, h_bulk, p, sa_final, ct_final, w_ih_final)
    sa_final = sa_final[]
    ct_final = ct_final[]
    w_ih_final = w_ih_final[]
    return sa_final, ct_final, w_ih_final
end

function gsw_frazil_properties_potential(sa_bulk, h_pot_bulk, p)
    sa_final = Ref{Float64}()
    ct_final = Ref{Float64}()
    w_ih_final = Ref{Float64}()
    ccall((:gsw_frazil_properties_potential, libgswteos), Cvoid, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa_bulk, h_pot_bulk, p, sa_final, ct_final, w_ih_final)
    sa_final = sa_final[]
    ct_final = ct_final[]
    w_ih_final = w_ih_final[]
    return sa_final, ct_final, w_ih_final
end

function gsw_frazil_properties_potential_poly(sa_bulk, h_pot_bulk, p)
    sa_final = Ref{Float64}()
    ct_final = Ref{Float64}()
    w_ih_final = Ref{Float64}()
    ccall((:gsw_frazil_properties_potential_poly, libgswteos), Cvoid, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa_bulk, h_pot_bulk, p, sa_final, ct_final, w_ih_final)
    sa_final = sa_final[]
    ct_final = ct_final[]
    w_ih_final = w_ih_final[]
    return sa_final, ct_final, w_ih_final
end

function gsw_frazil_ratios_adiabatic(sa, p, w_ih)
    dsa_dct_frazil = Ref{Float64}()
    dsa_dp_frazil = Ref{Float64}()
    dct_dp_frazil = Ref{Float64}()
    ccall((:gsw_frazil_ratios_adiabatic, libgswteos), Cvoid, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, p, w_ih, dsa_dct_frazil, dsa_dp_frazil, dct_dp_frazil)
    dsa_dct_frazil = dsa_dct_frazil[]
    dsa_dp_frazil = dsa_dp_frazil[]
    dct_dp_frazil = dct_dp_frazil[]
    return dsa_dct_frazil, dsa_dp_frazil, dct_dp_frazil
end

function gsw_frazil_ratios_adiabatic_poly(sa, p, w_ih)
    dsa_dct_frazil = Ref{Float64}()
    dsa_dp_frazil = Ref{Float64}()
    dct_dp_frazil = Ref{Float64}()
    ccall((:gsw_frazil_ratios_adiabatic_poly, libgswteos), Cvoid, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, p, w_ih, dsa_dct_frazil, dsa_dp_frazil, dct_dp_frazil)
    dsa_dct_frazil = dsa_dct_frazil[]
    dsa_dp_frazil = dsa_dp_frazil[]
    dct_dp_frazil = dct_dp_frazil[]
    return dsa_dct_frazil, dsa_dp_frazil, dct_dp_frazil
end

function gsw_ice_fraction_to_freeze_seawater(sa, ct, p, t_ih)
    sa_freeze = Ref{Float64}()
    ct_freeze = Ref{Float64}()
    w_ih = Ref{Float64}()
    ccall((:gsw_ice_fraction_to_freeze_seawater, libgswteos), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, t_ih, sa_freeze, ct_freeze, w_ih)
    sa_freeze = sa_freeze[]
    ct_freeze = ct_freeze[]
    w_ih = w_ih[]
    return sa_freeze, ct_freeze, w_ih
end

function gsw_ipv_vs_fnsquared_ratio(sa, ct, p, p_ref)
    nz = length(sa)
    ipv_vs_fnsquared_ratio = Array{Float64}(undef, nz - 1)
    p_mid = Array{Float64}(undef, nz - 1)
    ccall((:gsw_ipv_vs_fnsquared_ratio, libgswteos), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cint, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, p_ref, nz, ipv_vs_fnsquared_ratio, p_mid)
    return ipv_vs_fnsquared_ratio, p_mid
end

function gsw_linear_interp_sa_ct(sa, ct, p, p_i)
    np = length(p)
    npi = length(p_i)
    sa_i = Array{Float64}(undef, npi)
    ct_i = Array{Float64}(undef, npi)
    ccall((:gsw_linear_interp_sa_ct, libgswteos), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, np, p_i, npi, sa_i, ct_i)
    return sa_i, ct_i
end

function gsw_melting_ice_into_seawater(sa, ct, p, w_ih, t_ih)
    sa_final = Ref{Float64}()
    ct_final = Ref{Float64}()
    w_ih_final = Ref{Float64}()
    ccall((:gsw_melting_ice_into_seawater, libgswteos), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, w_ih, t_ih, sa_final, ct_final, w_ih_final)
    sa_final = sa_final[]
    ct_final = ct_final[]
    w_ih_final = w_ih_final[]
    return sa_final, ct_final, w_ih_final
end

function gsw_melting_seaice_into_seawater(sa, ct, p, w_seaice, sa_seaice, t_seaice)
    sa_final = Ref{Float64}()
    ct_final = Ref{Float64}()
    ccall((:gsw_melting_seaice_into_seawater, libgswteos), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, w_seaice, sa_seaice, t_seaice, sa_final, ct_final)
    sa_final = sa_final[]
    ct_final = ct_final[]
    return sa_final, ct_final
end

function gsw_pot_enthalpy_ice_freezing_first_derivatives(sa, p)
    pot_enthalpy_ice_freezing_sa = Ref{Float64}()
    pot_enthalpy_ice_freezing_p = Ref{Float64}()
    ccall((:gsw_pot_enthalpy_ice_freezing_first_derivatives, libgswteos), Cvoid, (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), sa, p, pot_enthalpy_ice_freezing_sa, pot_enthalpy_ice_freezing_p)
    pot_enthalpy_ice_freezing_sa = pot_enthalpy_ice_freezing_sa[]
    pot_enthalpy_ice_freezing_p = pot_enthalpy_ice_freezing_p[]
    return pot_enthalpy_ice_freezing_sa, pot_enthalpy_ice_freezing_p
end

function gsw_pot_enthalpy_ice_freezing_first_derivatives_poly(sa, p)
    pot_enthalpy_ice_freezing_sa = Ref{Float64}()
    pot_enthalpy_ice_freezing_p = Ref{Float64}()
    ccall((:gsw_pot_enthalpy_ice_freezing_first_derivatives_poly, libgswteos), Cvoid, (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), sa, p, pot_enthalpy_ice_freezing_sa, pot_enthalpy_ice_freezing_p)
    pot_enthalpy_ice_freezing_sa = pot_enthalpy_ice_freezing_sa[]
    pot_enthalpy_ice_freezing_p = pot_enthalpy_ice_freezing_p[]
    return pot_enthalpy_ice_freezing_sa, pot_enthalpy_ice_freezing_p
end

function gsw_pt_first_derivatives(sa, ct)
    pt_sa = Ref{Float64}()
    pt_ct = Ref{Float64}()
    ccall((:gsw_pt_first_derivatives, libgswteos), Cvoid, (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, pt_sa, pt_ct)
    pt_sa = pt_sa[]
    pt_ct = pt_ct[]
    return pt_sa, pt_ct
end

function gsw_pt_second_derivatives(sa, ct)
    pt_sa_sa = Ref{Float64}()
    pt_sa_ct = Ref{Float64}()
    pt_ct_ct = Ref{Float64}()
    ccall((:gsw_pt_second_derivatives, libgswteos), Cvoid, (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, pt_sa_sa, pt_sa_ct, pt_ct_ct)
    pt_sa_sa = pt_sa_sa[]
    pt_sa_ct = pt_sa_ct[]
    pt_ct_ct = pt_ct_ct[]
    return pt_sa_sa, pt_sa_ct, pt_ct_ct
end

function gsw_rho_alpha_beta(sa, ct, p)
    rho = Ref{Float64}()
    alpha = Ref{Float64}()
    beta = Ref{Float64}()
    ccall((:gsw_rho_alpha_beta, libgswteos), Cvoid, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, rho, alpha, beta)
    rho = rho[]
    alpha = alpha[]
    beta = beta[]
    return rho, alpha, beta
end

function gsw_rho_first_derivatives(sa, ct, p)
    drho_dsa = Ref{Float64}()
    drho_dct = Ref{Float64}()
    drho_dp = Ref{Float64}()
    ccall((:gsw_rho_first_derivatives, libgswteos), Cvoid, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, drho_dsa, drho_dct, drho_dp)
    drho_dsa = drho_dsa[]
    drho_dct = drho_dct[]
    drho_dp = drho_dp[]
    return drho_dsa, drho_dct, drho_dp
end

function gsw_rho_first_derivatives_wrt_enthalpy(sa, ct, p)
    rho_sa = Ref{Float64}()
    rho_h = Ref{Float64}()
    ccall((:gsw_rho_first_derivatives_wrt_enthalpy, libgswteos), Cvoid, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, rho_sa, rho_h)
    rho_sa = rho_sa[]
    rho_h = rho_h[]
    return rho_sa, rho_h
end

function gsw_rho_second_derivatives(sa, ct, p)
    rho_sa_sa = Ref{Float64}()
    rho_sa_ct = Ref{Float64}()
    rho_ct_ct = Ref{Float64}()
    rho_sa_p = Ref{Float64}()
    rho_ct_p = Ref{Float64}()
    ccall((:gsw_rho_second_derivatives, libgswteos), Cvoid, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, rho_sa_sa, rho_sa_ct, rho_ct_ct, rho_sa_p, rho_ct_p)
    rho_sa_sa = rho_sa_sa[]
    rho_sa_ct = rho_sa_ct[]
    rho_ct_ct = rho_ct_ct[]
    rho_sa_p = rho_sa_p[]
    rho_ct_p = rho_ct_p[]
    return rho_sa_sa, rho_sa_ct, rho_ct_ct, rho_sa_p, rho_ct_p
end

function gsw_rr68_interp_sa_ct(sa, ct, p, p_i)
    mp = length(p)
    mp_i = length(p_i)
    sa_i = Array{Float64}(undef, mp_i)
    ct_i = Array{Float64}(undef, mp_i)
    ccall((:gsw_rr68_interp_sa_ct, libgswteos), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, mp, p_i, mp_i, sa_i, ct_i)
    return sa_i, ct_i
end

function gsw_sa_freezing_estimate(p, saturation_fraction)
    t = Ref{Float64}()
    ct = Ref{Float64}()
    ccall((:gsw_sa_freezing_estimate, libgswteos), Cdouble, (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), p, saturation_fraction, ct, t)
    ct = ct[]
    t = t[]
    return ct, t
end

function gsw_seaice_fraction_to_freeze_seawater(sa, ct, p, sa_seaice, t_seaice)
    sa_freeze = Ref{Float64}()
    ct_freeze = Ref{Float64}()
    w_seaice = Ref{Float64}()
    ccall((:gsw_seaice_fraction_to_freeze_seawater, libgswteos), Cvoid, (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, sa_seaice, t_seaice, sa_freeze, ct_freeze, w_seaice)
    sa_freeze = sa_freeze[]
    ct_freeze = ct_freeze[]
    w_seaice = w_seaice[]
    return sa_freeze, ct_freeze, w_seaice
end

function gsw_specvol_alpha_beta(sa, ct, p)
    specvol = Ref{Float64}()
    alpha = Ref{Float64}()
    beta = Ref{Float64}()
    ccall((:gsw_specvol_alpha_beta, libgswteos), Cvoid, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, specvol, alpha, beta)
    specvol = specvol[]
    alpha = alpha[]
    beta = beta[]
    return specvol, alpha, beta
end

function gsw_specvol_first_derivatives(sa, ct, p)
    v_sa = Ref{Float64}()
    v_ct = Ref{Float64}()
    v_p = Ref{Float64}()
    ccall((:gsw_specvol_first_derivatives, libgswteos), Cvoid, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, v_sa, v_ct, v_p)
    v_sa = v_sa[]
    v_ct = v_ct[]
    v_p = v_p[]
    return v_sa, v_ct, v_p
end

function gsw_specvol_first_derivatives_wrt_enthalpy(sa, ct, p)
    v_sa = Ref{Float64}()
    v_h = Ref{Float64}()
    ccall((:gsw_specvol_first_derivatives_wrt_enthalpy, libgswteos), Cvoid, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, v_sa, v_h)
    v_sa = v_sa[]
    v_h = v_h[]
    return v_sa, v_h
end

function gsw_specvol_second_derivatives(sa, ct, p)
    v_sa_sa = Ref{Float64}()
    v_sa_ct = Ref{Float64}()
    v_ct_ct = Ref{Float64}()
    v_sa_p = Ref{Float64}()
    v_ct_p = Ref{Float64}()
    ccall((:gsw_specvol_second_derivatives, libgswteos), Cvoid, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, v_sa_sa, v_sa_ct, v_ct_ct, v_sa_p, v_ct_p)
    v_sa_sa = v_sa_sa[]
    v_sa_ct = v_sa_ct[]
    v_ct_ct = v_ct_ct[]
    v_sa_p = v_sa_p[]
    v_ct_p = v_ct_p[]
    return v_sa_sa, v_sa_ct, v_ct_ct, v_sa_p, v_ct_p
end

function gsw_specvol_second_derivatives_wrt_enthalpy(sa, ct, p)
    v_sa_sa = Ref{Float64}()
    v_sa_h = Ref{Float64}()
    v_h_h = Ref{Float64}()
    ccall((:gsw_specvol_second_derivatives_wrt_enthalpy, libgswteos), Cvoid, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, v_sa_sa, v_sa_h, v_h_h)
    v_sa_sa = v_sa_sa[]
    v_sa_h = v_sa_h[]
    v_h_h = v_h_h[]
    return v_sa_sa, v_sa_h, v_h_h
end

function gsw_t_freezing_first_derivatives(sa, p, saturation_fraction)
    tfreezing_sa = Ref{Float64}()
    tfreezing_p = Ref{Float64}()
    ccall((:gsw_t_freezing_first_derivatives, libgswteos), Cvoid, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), sa, p, saturation_fraction, tfreezing_sa, tfreezing_p)
    tfreezing_sa = tfreezing_sa[]
    tfreezing_p = tfreezing_p[]
    return tfreezing_sa, tfreezing_p
end

function gsw_t_freezing_first_derivatives_poly(sa, p, saturation_fraction)
    tfreezing_sa = Ref{Float64}()
    tfreezing_p = Ref{Float64}()
    ccall((:gsw_t_freezing_first_derivatives_poly, libgswteos), Cvoid, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), sa, p, saturation_fraction, tfreezing_sa, tfreezing_p)
    tfreezing_sa = tfreezing_sa[]
    tfreezing_p = tfreezing_p[]
    return tfreezing_sa, tfreezing_p
end


function gsw_util_interp1q_int(x, iy)
    nx = length(x)
    nxi = length(iy)
    x_i = Array{Float64}(undef, nxi)
    y_i = Array{Float64}(undef, nxi)
    ccall((:gsw_util_interp1q_int, libgswteos), Ptr{Cdouble}, (Cint, Ptr{Cdouble}, Ptr{Cint}, Cint, Ptr{Cdouble}, Ptr{Cdouble}), nx, x, iy, nxi, x_i, y_i)
    return x_i, y_i
end

function gsw_util_linear_interp(x, y)
    nx = length(x)
    ny = length(y)
    nxi = length(x_i)
    x_i = Array{Float64}(undef, nxi)
    y_i = Array{Float64}(undef, nxi)
    ccall((:gsw_util_linear_interp, libgswteos), Ptr{Cdouble}, (Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}), nx, x, ny, y, nxi, x_i, y_i)
    return x_i, y_i
end

function gsw_util_sort_real(rarray)
    nx = length(rarray)
    iarray = Array{Int32}(undef, nx)
    ccall((:gsw_util_sort_real, libgswteos), Cvoid, (Ptr{Cdouble}, Cint, Ptr{Cint}), rarray, nx, iarray)
    return iarray
end

function gsw_util_xinterp1(x, y)
    n = length(x)
    x0 = Ref{Float64}(x0)
    ccall((:gsw_util_xinterp1, libgswteos), Cdouble, (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble), x, y, n, x0)
    x0 = x0[]
    return x0
end

function gsw_util_pchip_interp(x, y)
    n = length(x)
    ni = length(xi)
    xi = Array{Float64}(undef, ni)
    yi = Array{Float64}(undef, ni)
    ccall((:gsw_util_pchip_interp, libgswteos), Cint, (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint), x, y, n, xi, yi, ni)
    return xi, yi
end
