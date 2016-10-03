module FunctionZeros

using Roots

export besselj_zero

# Asymptotic formula for zeros of Bessel J function of order nu
"""
    besselj_asymptotic_zero(nu,n)

Asymptotic formula for the `n`th zero fo the the Bessel J function of order `nu`.
"""
function besselj_asymptotic_zero(nu,n)
    return pi * (n-1 + nu/2 + 3//4)
end

# Find nth zero of Bessel J. Use the asymptotic values as
# starting values. These find the correct zeros even
# for n = 1,2,...

"""
    besselj_zero(nu,n)

`n`th zero of the Bessel J function of order `nu`,
for `n` = `1,2,...`.
"""
function besselj_zero(nu,n)
    z = besselj_asymptotic_zero(nu,n)
    bf = (x) -> besselj(nu,x)
    return fzero(bf,z,ftol=1e-15)
end

end # module FunctionZeros
