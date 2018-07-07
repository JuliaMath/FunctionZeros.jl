module FunctionZeros
import SpecialFunctions
import Roots

export besselj_zero

# Asymptotic formula for zeros of Bessel J function of order nu
"""
    besselj_asymptotic_zero(nu,n)

Asymptotic formula for the `n`th zero fo the the Bessel J function of order `nu`.
`besselj_asymptotic_zero` is vectorized.
"""
besselj_asymptotic_zero(nu, n) = pi * (n - 1 + nu / 2 + 3//4)

# Use the asymptotic values as starting values.
# These find the correct zeros even for n = 1,2,...
# Order 0 is 6 times slower and 50-100 times less accurate,
# with other parameters constant that higher orders.
"""
    besselj_zero(nu,n)

`n`th zero of the Bessel J function of order `nu`,
for `n` = `1,2,...`.
"""
besselj_zero(nu,n) = Roots.fzero((x) -> SpecialFunctions.besselj(nu, x),
                                 besselj_asymptotic_zero(nu, n); order=1)

end # module FunctionZeros
