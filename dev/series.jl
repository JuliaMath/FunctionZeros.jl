# This is an asymptotic form with more terms than
# the series used in ../src/. It gives no
# advantage as the initial point for root finding
# since the root finding takes far more time.
#
# This could be used for large enough n without
# root finding, especially if the user specifies
# a desired precision.
function Jzeros(m::Real,n::Integer)
    beta = MathConstants.pi * (n + m / 2 - 1//4)
    delta = 8 * beta
    mu = 4 * m^2
    t1 = 1
    t2 = 4 * (7 * mu - 31) / (3 * delta^2)
    t3 = 32 * (84 * mu^2 - 982 * mu + 3779) / (15 * delta^4)
    t4 = 64 * (6949 * mu^3 - 153855 * mu^2 + 1585743 * mu - 6277237) /
        (105 * delta^6)
    return beta - (mu - 1) / delta * (t1 + t2 + t3 + t4)
end
