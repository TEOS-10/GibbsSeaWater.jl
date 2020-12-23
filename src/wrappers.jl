function gsw_nsquared(sa::AbstractArray{T}, ct::AbstractArray{T}, p::AbstractArray{T}, lat::AbstractArray{T}) where {T <: Real}
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