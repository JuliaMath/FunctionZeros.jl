# FunctionZeros
*Zeros of the Bessel J and Y functions and their derivatives*

[![Build Status](https://github.com/JuliaMath/FunctionZeros.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/JuliaMath/FunctionZeros.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Codecov](https://codecov.io/gh/JuliaMath/ILog2.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaMath/FunctionZeros.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![JET QA](https://img.shields.io/badge/JET.jl-%E2%9C%88%EF%B8%8F-%23aa4444)](https://github.com/aviatesk/JET.jl)

This package provides functions to compute the zeros of the J and Y functions,
and the zeros of their derivatives, where J and Y are Bessel functions of the first and second kind, respectively. 

For all functions described below, the order `nu::Real` is a finite number and `n::Integer` is a positive integer.
When `nu isa AbstractFloat`, the returned value has the same type as `nu`. When `nu isa Integer`, the usual
promotion rules apply, so that for most builtin integer types the output type will be `Float64`. However, 
when `nu isa BigInt` the output type will be `BigFloat`.

When the output type is `Float64`, the exported functions (`besselj_zero`, 
`bessely_zero`, `besselj_deriv_zero`, and `bessely_deriv_zero`) will use lookup tables to rapidly
return function zeros if the order `nu` is one of the first few values of `0, 1, ...` and the enumerator
`n` is one of the first values of `1, 2, 3, ...`.  See the individual function docstrings for the actual
extents of the lookup tables.

### Exported Functions

#### besselj_zero(nu, n)

```julia
besselj_zero(nu, n)
```

Return the `n`th zero of the Bessel J function of order `nu`. 

#### bessely_zero(nu, n)

```julia
bessely_zero(nu, n)
```

Return the `n`th zero of the Bessel Y function of order `nu`.

#### besselj_deriv_zero(nu, n)

```julia
besselj_deriv_zero(nu, n)
```

Return the `n`th nonvanishing zero of the derivative of the Bessel J
function of order `nu`.

#### bessely_deriv_zero(nu, n)

```julia
bessely_deriv_zero(nu, n)
```

Return the `n`th zero of the derivative of the Bessel Y function of order `nu`.

### Non-exported But Useful Functions

#### FunctionZeros.besselj_zero_asymptotic(nu, n)

```julia
FunctionZeros.besselj_zero_asymptotic(nu, n)
```

Asymptotic formula for the `n`th zero for the Bessel J function of order `nu`.


#### FunctionZeros.bessely_zero_asymptotic(nu, n)

```julia
FunctionZeros.bessely_zero_asymptotic(nu, n)
```

Asymptotic formula for the `n`th zero for the Bessel Y function of order `nu`.


#### FunctionZeros.besselj_deriv_zero_asymptotic(nu, n)

```julia
FunctionZeros.besselj_deriv_zero_asymptotic(nu, n)
```

Asymptotic formula for the `n`th nonvanishing zero of the derivative of the 
Bessel J function of order `nu`.


#### FunctionZeros.bessely_deriv_zero_asymptotic(nu, n)

```julia
FunctionZeros.bessely_deriv_zero_asymptotic(nu, n)
```

Asymptotic formula for the `n`th zero of the derivative of the Bessel Y function of order `nu`.


