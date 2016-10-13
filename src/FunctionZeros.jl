VERSION >= v"0.4.0-dev+6521" && __precompile__()

module FunctionZeros

using Roots

export besselj_zero

# Asymptotic formula for zeros of Bessel J function of order nu
"""
    besselj_asymptotic_zero(nu,n)

Asymptotic formula for the `n`th zero fo the the Bessel J function of order `nu`.
`besselj_asymptotic_zero` is vectorized.
"""
besselj_asymptotic_zero(nu,n) = return pi * (n-1 + nu/2 + 3//4)

# Find nth zero of Bessel J. Use the asymptotic values as
# starting values. These find the correct zeros even
# for n = 1,2,...

"""
    besselj_zero(nu,n)

`n`th zero of the Bessel J function of order `nu`,
for `n` = `1,2,...`. `besselj_zero` is vectorized.
"""
function besselj_zero(nu,n)
    z = besselj_asymptotic_zero(nu,n)
    bf = (x) -> besselj(nu,x)
#    return fzero(bf,z,ftol=1e-15)
    return fzero(bf,z)    
end

Base.@vectorize_2arg Number besselj_zero
Base.@vectorize_2arg Number besselj_asymptotic_zero

end # module FunctionZeros
