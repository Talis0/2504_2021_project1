#############################################################################
#############################################################################
#
# This file implements polynomial GCD 
#                                                                               
#############################################################################
#############################################################################

"""
The extended euclid algorithm for polynomials modulo prime.
"""
function extended_euclid_alg(a::Polynomial, b::Polynomial, prime::Int)
    old_r, r = mod(a,prime), mod(b,prime)
    old_s, s = one(Polynomial), zero(Polynomial)
    old_t, t = zero(Polynomial), one(Polynomial)

    while !iszero(mod(r,prime))

        q = divide(old_r, r)(prime) |> first

        old_r, r = r, mod(old_r - q*r, prime)

        old_s, s = s, mod(old_s - q*s, prime)

        old_t, t = t, mod(old_t - q*t, prime)

    end
    g, s, t = old_r, old_s, old_t
    @assert mod(s*a + t*b - g, prime) == 0
    return g, s, t  
end

function extended_euclid_alg(a::PolynomialModP, b::PolynomialModP)
    prime = a.mod
    old_r, r = a, b
    old_s, s = one(PolynomialModP,prime), zero(PolynomialModP,prime)
    old_t, t = zero(PolynomialModP,prime), one(PolynomialModP,prime)

    while !iszero(r)

        q = divide(old_r, r) |> first

        old_r, r = r, old_r - q*r

        old_s, s = s, old_s - q*s

        old_t, t = t, old_t - q*t

    end
    g, s, t = old_r, old_s, old_t
    @assert s*a + t*b - g == 0
    return g, s, t  
end

"""
The GCD of two polynomials modulo prime.
"""
#= gcd(a::Polynomial, b::Polynomial, prime::Int) = extended_euclid_alg(a,b,prime) |> first =#
gcd(a::PolynomialModP, b::PolynomialModP) = extended_euclid_alg(a,b) |> first