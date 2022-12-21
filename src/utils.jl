Base.@irrational Ωe 0.00007292115 big(0.00007292115)
Ωe, const Omegae = Ωe
gsw_f(latitude::T) where {T<:Real} = T(2) * Ωe * sind(latitude)

# constants from gsw_internal_const.h

# cp0  =  The "specific heat" for use    [ J/(kg K) ]\
#           with Conservative Temperature   \

Base.@irrational gsw_cp0 3991.86795711963 big(3991.86795711963)
# gsw_cp0, const gsw_cp0

#  T0  =  the Celsius zero point.                                     [ K ]\
Base.@irrational gsw_t0 273.15 big(273.15)

# P0  =  Absolute Pressure of one standard atmosphere.              [ Pa ]\
Base.@irrational gsw_p0 101325.0 big(101325.0)

#  SSO  =  Standard Ocean Reference Salinity.                      [ g/kg ]\
Base.@irrational gsw_sso 35.16504 big(35.16504)
Base.@irrational gsw_sqrtsso 5.930011804372737 big(5.930011804372737)

#  uPS  =  unit conversion factor for salinities                   [ g/kg ]\
Base.@irrational gsw_ups 1.00471542857142850380 big(1.00471542857142850380) # gsw_sso / 35.0

#  sfac  =  1/(40*gsw_ups)\
Base.@irrational gsw_sfac 0.0248826675584615 big(0.0248826675584615)

#  deltaS = 24, offset = deltaS*gsw_sfac\
Base.@irrational offset 0.5971840214030754 big(0.5971840214030754)

#  C3515  =  Conductivity at (SP=35, t_68=15, p=0)                [ mS/cm ]\
Base.@irrational gsw_c3515 42.9140 big(42.9140)

#  SonCl  =  SP to Chlorinity ratio                           [ (g/kg)^-1 ]\
Base.@irrational gsw_soncl 1.80655 big(1.80655)


#  valence_factor  =  valence factor of sea salt of Reference Composition\
#                                                             [ unitless ]\
Base.@irrational gsw_valence_factor 1.2452898 big(1.2452898)

#  atomic_weight = mole-weighted atomic weight of sea salt of Reference \
#                  Composition                                    [ g/mol ]\
Base.@irrational gsw_atomic_weight 31.4038218 big(31.4038218)
