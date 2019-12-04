using GSW
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
    as = 35.0; ct = 10.0; pr = 10.0
    println(GSW.gsw_rho(as, ct, pr))
end

# write your own tests here
#@test 1 == 2
