# FunctionZeros

This module provides a function to compute the zeros of the Bessel J function.
No tests are written.

#### besselj_zero

```julia
besselj_zero(nu,n)
```

Return the `n`th zero of the the Bessel J function of order `nu`. The returned
type has the same type as `nu`.


#### Notes

There are a few modules like this one implementing some special functions. They should
probably be collected in some organized way.
