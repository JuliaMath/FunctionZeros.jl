# FunctionZeros
*Zeros of the Bessel J function*

Linux, OSX: [![Build Status](https://travis-ci.org/jlapeyre/FunctionZeros.jl.svg)](https://travis-ci.org/jlapeyre/FunctionZeros.jl)
&nbsp;
Windows: [![Build Status](https://ci.appveyor.com/api/projects/status/github/jlapeyre/FunctionZeros.jl?branch=master&svg=true)](https://ci.appveyor.com/project/jlapeyre/functionzeros-jl)
&nbsp; &nbsp; &nbsp;
[![Coverage Status](https://coveralls.io/repos/jlapeyre/FunctionZeros.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/jlapeyre/FunctionZeros.jl?branch=master)
[![codecov.io](http://codecov.io/github/jlapeyre/FunctionZeros.jl/coverage.svg?branch=master)](http://codecov.io/github/jlapeyre/FunctionZeros.jl?branch=master)

This module provides a function to compute the zeros of the Bessel J function.

#### besselj_zero(nu,n)

```julia
besselj_zero(nu,n)
```

Return the `n`th zero of the the Bessel J function of order `nu`. The returned
type has the same type as `nu`.

#### besselj_asymptotic_zero(nu,n)

Use asymptotic form expression for the zeros. This function is not exported.

#### Notes

There are a few modules like this one implementing some special functions. They should
probably be collected in some organized way.
