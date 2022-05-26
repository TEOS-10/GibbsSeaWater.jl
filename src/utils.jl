gsw_f(latitude::T) where {T<:Real} = T(2) * T(7.292115e-5) * sin(latitude * pi / T(180.0))
