module FunctionZeros
import SpecialFunctions
import Roots

export besselj_zero, bessely_zero

"""
    besselj_zero_asymptotic_four_terms(nu, n, kind=1)

Asymptotic formula for the `n`th zero of the the Bessel J function of order `nu`.
`kind == 1 (2)` for Bessel function of the first (second) kind, J (Y).
"""
function besselj_zero_asymptotic_four_terms(nu::Real, n::Integer, kind=1)
    if kind == 1
        beta = MathConstants.pi * (n + nu / 2 - 1//4)
    else # kind == 2
        beta = MathConstants.pi * (n + nu / 2 - 3//4)
    end
    delta = 8 * beta
    mu = 4 * nu^2
    mup2 = mu * mu
    mup3 = mup2 * mu
    deltap2 = delta * delta
    deltap3 = deltap2 * delta
    deltap4 = deltap2 * deltap2
    deltap6 = deltap3 * deltap3
    t1 = 1
    t2 = 4 * (7 * mu - 31) / (3 * deltap2)
    t3 = 32 * (83 * mup2 - 982 * mu + 3779) / (15 * deltap4)
    t4 = 64 * (6949 * mup3 - 153855 * mup2 + 1585743 * mu - 6277237) /
        (105 * deltap6)
    zero_asymp = beta - ((mu - 1) / delta) * (t1 + t2 + t3 + t4)
    return zero_asymp
end


# Asymptotic formula for zeros of Bessel J function of order nu
"""
    besselj_zero_asymptotic(nu, n)

Asymptotic formula for the `n`th zero of the the Bessel J function of order `nu`.
"""
besselj_zero_asymptotic(nu, n) = pi * (n - 1 + nu / 2 + 3//4)

"""
    bessely_zero_asymptotic(nu, n)

Asymptotic formula for the `n`th zero of the the Bessel Y function of order `nu`.
"""
bessely_zero_asymptotic(nu, n) = pi * (n + nu / 2 - 3//4)

# Use the asymptotic values as starting values.
# These find the correct zeros even for n = 1,2,...
# Order 0 is 6 times slower and 50-100 times less accurate
# than higher orders, with other parameters constant.
"""
    besselj_zero(nu, n; order=2)

`n`th zero of the Bessel J function of order `nu`,
for `n` = `1,2,...`.

`order` is passed to the function `Roots.fzero`.
"""
besselj_zero(nu, n; order=2) = Roots.fzero((x) -> SpecialFunctions.besselj(nu, x),
                                           besselj_zero_asymptotic_four_terms(nu, n); order=order)

"""
    bessely_zero(nu, n; order=2)

`n`th zero of the Bessel Y function of order `nu`,
for `n` = `1,2,...`.

`order` is passed to the function `Roots.fzero`.
"""
bessely_zero(nu, n; order=2) = Roots.fzero((x) -> SpecialFunctions.bessely(nu, x),
                                           besselj_zero_asymptotic_four_terms(nu, n, 2); order=order)

end # module FunctionZeros
