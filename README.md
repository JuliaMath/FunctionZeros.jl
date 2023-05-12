# FunctionZeros
*Zeros of the Bessel J and Y functions*

[![Build Status](https://github.com/JuliaMath/FunctionZeros.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/JuliaMath/FunctionZeros.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage Status](https://coveralls.io/repos/github/JuliaMath/FunctionZeros.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaMath/FunctionZeros.jl?branch=master)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![JET QA](https://img.shields.io/badge/JET.jl-%E2%9C%88%EF%B8%8F-%23aa4444)](https://github.com/aviatesk/JET.jl)

This module provides a function to compute the zeros of the Bessel J and K functions,
that is Bessel functions of the first and second kind.

#### besselj_zero(nu, n)

```julia
besselj_zero(nu, n)
```

Return the `n`th zero of the the Bessel J function of order `nu`. The returned
type has the same type as `nu`.

#### FunctionZeros.besselj_zero_asymptotic(nu, n)

Asymptotic formula for the `n`th zero fo the the Bessel J function of order `nu`.

```julia
besselj_zero_asymptotic(nu, n)
```


#### bessely_zero(nu, n)

```julia
bessely_zero(nu, n)
```

Return the `n`th zero of the the Bessel Y function of order `nu`. The returned
type has the same type as `nu`.

#### FunctionZeros.bessely_zero_asymptotic(nu, n)

Asymptotic formula for the `n`th zero fo the the Bessel Y function of order `nu`.

```julia
bessely_zero_asymptotic(nu, n)
```
