
"""
    FunctionZeros
This module provides functions to compute the zeros of the J and Y functions,
and the zeros of their derivatives, where J and Y are Bessel functions of the first and second kind, respectively.
"""
module FunctionZeros
import SpecialFunctions
import Roots

export besselj_zero, bessely_zero, besselj_deriv_zero, bessely_deriv_zero

# Set max order and index of zeros to precompute and tabulate
const nupre_max = 1
const npre_max = 500

# Strings used in multiple function docstrings:
const speeddocstr = """For greater speed, table lookup is used for `Float64` outputs when
                       `nu ∈ 0:$nupre_max` and `n ∈ 1:$(npre_max)`."""
const argstr = """## Arguments
                  - `nu::Real`: The order of the Bessel function.
                  - `n::Integer`: Enumerates the zero to be found.

                  ## Return Value
                  The return value type is determined by `nu`.
                  When `nu isa AbstractFloat`, the returned value has the same type as `nu`.
                  When `nu isa Integer`, the usual promotion rules apply.
               """

"""
    besselj_zero_asymptotic(nu, n)

Asymptotic formula for the `n`th zero of the the Bessel function of the first kind J of order `nu`.

$argstr
"""
besselj_zero_asymptotic(nu, n) = bessel_zero_asymptotic(nu, n, 1)

"""
    bessely_zero_asymptotic(nu, n)

Asymptotic formula for the `n`th zero of the the Bessel function of the second kind Y of order `nu`.

$argstr
"""
bessely_zero_asymptotic(nu, n) = bessel_zero_asymptotic(nu, n, 2)


"""
    bessel_zero_asymptotic(nu, n, kind=1)

Asymptotic formula for the `n`th zero of the the Bessel J (Y) function of order `nu`.
`kind == 1 (2)` for Bessel function of the first (second) kind, J (Y).
"""
function bessel_zero_asymptotic(nu_in::Real, n::Integer, kind=1)
    nu = abs(nu_in)
    if nu > 25 && n ≤ length(airy_zeros().ai)
        return bessel_zero_largenu_asymptotic(nu, n, kind, false)
    end
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

bessel_zero_asymptotic(nu::Integer, n::Integer, kind=1) = bessel_zero_asymptotic(float(nu), n, kind) # fix #27

# Use the asymptotic values as starting values.
# These find the correct zeros even for n = 1,2,...
# Order 0 is 6 times slower and 50-100 times less accurate
# than higher orders, with other parameters constant.
"""
    _besselj_zero(nu, n)

`n`th zero of the Bessel J function of order `nu`,
for `n` = `1,2,...`.

$argstr
"""
function _besselj_zero(nu::Real, n::Integer)
    return Roots.find_zero(bessel_zero_asymptotic(nu, n, 1)) do x
        return SpecialFunctions.besselj(nu, x)
    end
end

# Float64 tabulation of selected besselj_zero values
const jzero_pre = [_besselj_zero(nu, n) for nu in 0:nupre_max, n in 1:npre_max]

"""
    besselj_zero(nu, n)

Return the `n`th zero of the Bessel J function of order `nu`,
for `n` = `1,2,...`.

$argstr
$speeddocstr
"""
besselj_zero(nu::Real, n::Integer) = _besselj_zero(nu, n)

besselj_zero(nu::BigInt, n::Integer) = _besselj_zero(nu, n)

function besselj_zero(nu::Union{Integer,Float64}, n::Integer)
    if nu in 0:nupre_max && n in 1:npre_max
        return jzero_pre[Int(nu) + 1, Int(n)]
    else
        return _besselj_zero(nu, n)
    end
end


"""
    _bessely_zero(nu, n)

`n`th zero of the Bessel Y function of order `nu`,
for `n` = `1,2,...`.

$argstr
"""
function _bessely_zero(nu, n)
    if isone(n) && abs(nu) < 0.1587 # See Issue 21
        return Roots.find_zero((nu + besselj_zero(nu, n)) / 2) do x
            SpecialFunctions.bessely(nu, x)
        end
    else
        return Roots.find_zero(bessel_zero_asymptotic(nu, n, 2)) do x
            SpecialFunctions.bessely(nu, x)
        end
    end
end

# tabulation of selected bessely_zero values in Float64
const yzero_pre = [_bessely_zero(nu, n) for nu in 0:nupre_max, n in 1:npre_max]

"""
    bessely_zero(nu, n)

Return the `n`th zero of the Bessel Y function of order `nu`,
for `n` = `1,2,...`.

$argstr
$speeddocstr
"""
bessely_zero(nu::Real, n::Integer) = _bessely_zero(nu, n)

bessely_zero(nu::BigInt, n::Integer) = _bessely_zero(nu, n)

function bessely_zero(nu::Union{Integer,Float64}, n::Integer)
    if nu in 0:nupre_max && n in 1:npre_max
        return yzero_pre[Int(nu) + 1, Int(n)]
    else
        return _bessely_zero(nu, n)
    end
end

"""
    besselj_deriv_zero_asymptotic(nu, n)

Asymptotic formula for the `n`th zero of the derivative of the Bessel function of the first kind J of order `nu`.

$argstr
"""
besselj_deriv_zero_asymptotic(nu, n) = bessel_deriv_zero_asymptotic(nu, n, 1)

"""
    bessely_deriv_zero_asymptotic(nu, n)

Asymptotic formula for the `n`th zero of the derivative of the Bessel function of the second kind Y of order `nu`.

$argstr
"""
bessely_deriv_zero_asymptotic(nu, n) = bessel_deriv_zero_asymptotic(nu, n, 2)


"""
    bessel_deriv_zero_asymptotic(nu, n, kind=1)

Asymptotic formula for the `n`th zero of the the derivative of Bessel J (Y) function
of order `nu`. `kind == 1 (2)` for Bessel function of the first (second) kind, J (Y).
"""
function bessel_deriv_zero_asymptotic(nu_in::Real, n::Integer, kind=1)
    nu = abs(nu_in)
    if nu > 25 && n ≤ length(airy_zeros().ai)
        return bessel_zero_largenu_asymptotic(nu, n, kind, true)
    end

    # Reference: https://dlmf.nist.gov/10.21.E20
    if kind == 1
        beta = MathConstants.pi * (n + nu / 2 - 3//4)
    else # kind == 2
        beta = MathConstants.pi * (n + nu / 2 - 1//4)
    end
    delta = 8 * beta
    mu = 4 * nu^2
    mup2 = mu * mu
    mup3 = mup2 * mu
    mup4 = mup3 * mu
    deltap2 = delta * delta
    deltap3 = deltap2 * delta
    deltap5 = deltap3 * deltap2
    deltap7 = deltap5 * deltap2
    t1 = (mu + 3) / delta
    t2 = 4 * (7 * mup2 + 82 * mu - 9) / (3 * deltap3)
    t3 = 32 * (83 * mup3 + 2075 * mup2 - 3039 * mu + 3537) / (15 * deltap5)
    t4 = 64 * (6949 * mup4 + 296492 * mup3 - 1248002 * mup2 + 7414380 * mu - 5853627) /
        (105 * deltap7)
    zero_asymp = beta - (t1 + t2 + t3 + t4)
    return zero_asymp
end

bessel_deriv_zero_asymptotic(nu::Integer, n::Integer, kind=1) = bessel_deriv_zero_asymptotic(float(nu), n, kind) # fix #27

"""
    _besselj_deriv_zero(nu, n)

Return the `n`th nonvanishing zero of the derivative of Bessel J function of order `nu`,
for `n` = `1,2,...`.

$argstr
"""
function _besselj_deriv_zero(nu::Real, n::Integer)
    # Ref: https://dlmf.nist.gov/10.6.E1
    iszero(nu) && (n += 1) # Skip the zero occuring at zero
    return Roots.find_zero(bessel_deriv_zero_asymptotic(nu, n, 1)) do x
        SpecialFunctions.besselj(nu - 1, x) - SpecialFunctions.besselj(nu + 1, x)
    end
end

# tabulation of selected besselj_deriv_zero values in Float64
const jderivzero_pre = [_besselj_deriv_zero(nu, n) for nu in 0:nupre_max, n in 1:npre_max]

"""
    besselj_deriv_zero(nu, n)

Return the `n`th nonvanishing zero of the derivative of the Bessel J function of order `nu`,
for `n` = `1,2,...`.

$argstr
$speeddocstr
"""
besselj_deriv_zero(nu::Real, n::Integer) = _besselj_deriv_zero(nu, n)

besselj_deriv_zero(nu::BigInt, n::Integer) = _besselj_deriv_zero(nu, n)

function besselj_deriv_zero(nu::Union{Integer,Float64}, n::Integer)
    if nu in 0:nupre_max && n in 1:npre_max
        return jderivzero_pre[Int(nu) + 1, Int(n)]
    else
        return _besselj_deriv_zero(nu, n)
    end
end

"""
    _bessely_deriv_zero(nu, n)

Return the `n`th zero of the derivative of the Bessel Y function of order `nu`,
for `n` = `1,2,...`.

$argstr
"""
function _bessely_deriv_zero(nu::Real, n::Integer)
    # Ref: https://dlmf.nist.gov/10.6.E1
    return Roots.find_zero(bessel_deriv_zero_asymptotic(nu, n, 2)) do x
        SpecialFunctions.bessely(nu - 1, x) - SpecialFunctions.bessely(nu + 1, x)
    end
end

# tabulation of selected bessely_deriv_zero values in Float64
const yderivzero_pre = [_bessely_deriv_zero(nu, n) for nu in 0:nupre_max, n in 1:npre_max]

"""
    bessely_deriv_zero(nu, n)

Return the `n`th zero of the derivative of the Bessel Y function of order `nu`,
for `n` = `1,2,...`.

$argstr
$speeddocstr
"""
bessely_deriv_zero(nu::Real, n::Integer) = _bessely_deriv_zero(nu, n)

bessely_deriv_zero(nu::BigInt, n::Integer) = _bessely_deriv_zero(nu, n)

function bessely_deriv_zero(nu::Union{Integer,Float64}, n::Integer)
    if nu in 0:nupre_max && n in 1:npre_max
        return yderivzero_pre[Int(nu) + 1, Int(n)]
    else
        return _bessely_deriv_zero(nu, n)
    end
end

"""
    airy_zeros()

Return the first few negative zeros of the functions `airyai`, `airybi`, `airyaiprime`, and `airybiprime`

## Return Value
- `(; ai, bi, aiprime, biprime)`: A named tuple containing the first few negative zeros of the functions
  `airyai`, `airybi`, `airyaiprime`, and `airybiprime`, respectively, as defined in the `SpecialFunctions`
   package. Each field in the named tuple consists of a tuple of `n` increasingly negative values. Here `n`
   is a small integer, say 2 or 3. 
"""
@inline function airy_zeros()
    ai = (-2.338107410459767, -4.087949444130973, -5.520559828095556)
    bi = (-1.173713222709127, -3.2710933028363516, -4.8307378416620095)
    aiprime = (-1.0187929716474717, -3.248197582179841, -4.820099211178737)
    biprime = (-2.2944396826141227, -4.073155089071816, -5.5123957296635835)
    return (; ai, bi, aiprime, biprime)
end

"""
    bessel_zero_largenu_asymptotic(nu::Real, m::Integer, kind::Integer, deriv::Bool)

Compute an asymptotic approximation for a zero of the J or Y Bessel function (or their derivative) suitable for large order.

## Input Arguments
- `nu`: The (nonnegative) order of the Bessel function.
- `m`: A small positive integer enumerating which zero is sought. Can not exceed `length(airy_zeros().ai)`.
  The asymptotic formulas used herein weaken as `m` increases.
- `kind`: Kind of Bessel function whose zero is sought. Either 1 for the J Bessel function or 2 for the Y Bessel function.
- `deriv`: If false, returns the zero of the function.  If true, returns the zero of the function's derivative.

## Return Value
`z`: A `Float64` approximation of the desired zero.

## Reference
- https://dlmf.nist.gov/10.21.vii       
"""
function bessel_zero_largenu_asymptotic(nu::Real, m::Integer, kind::Integer, deriv::Bool)
    abfactor = inv(cbrt(2.0))

    if kind == 1
        alpha = -abfactor * (deriv ? airy_zeros().aiprime[m] : airy_zeros().ai[m])
    elseif kind == 2
        alpha = -abfactor * (deriv ? airy_zeros().biprime[m] : airy_zeros().bi[m])
    else
        throw(ArgumentError("kind must be 1 or 2 but is $kind"))
    end

    alphap2 = alpha * alpha
    alphap3 = alpha * alphap2
    alphap4 = alpha * alphap3
    alphap5 = alpha * alphap4

    alpha0 = 1.0
    alpha1 = alpha

    if deriv
        alpha2 = 3 * alphap2 / 10 - inv(10 * alpha)
        alpha3 = -(alphap3 / 350 + inv(25) + inv(200 * alphap3))
        alpha4 = -479 * alphap4 / 630_000 + 509 * alpha / 31_500 + inv(1500 * alphap2) - inv(2_000 * alphap5)
        alpha5 = 0.0
    else
        alpha2 = 3 * alphap2 / 10
        alpha3 = -alphap3 / 350 + inv(70)
        alpha4 = -(479 * alphap4 / 63_000 + alpha / 3150)
        alpha5 = 20_231 * alphap5 / 8_085_000 - 551 * alphap2 / 161_700
    end

    x = inv(cbrt(nu)^2)
    xpk = 1.0
    zsum = lastterm = alpha0
    for alphak in (alpha1, alpha2, alpha3, alpha4, alpha5)
        xpk *= x
        nextterm = alphak * xpk
        abs(nextterm) > abs(lastterm) && break # Asymptotic series starting to diverge
        zsum += nextterm
        lastterm = nextterm
    end

    zsum *= nu
    return zsum
end

end # module FunctionZeros
