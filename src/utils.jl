Base.@irrational Ωe 0.00007292115 big(0.00007292115)
Ωe, const Omegae = Ωe
gsw_f(latitude::T) where {T<:Real} = T(2) * Ωe * sind(latitude)
