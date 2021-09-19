#############################################################################
#############################################################################
#
# This file implements polynomial multiplication 
#                                                                               
#############################################################################
#############################################################################

"""
Multiply two polynomials.
"""
function *(p1::Polynomial, p2::Polynomial)::Polynomial
    p_out = Polynomial()
    for t in p1
        p_out = p_out + (t * p2)
    end
    return p_out
end

function *(p1::PolynomialModP, p2::PolynomialModP)::PolynomialModP
    @assert p1.mod == p2.mod
    p = p1.mod
    p_out = PolynomialModP(p)
    for t in p1
        p_out = p_out + (t * p2)
    end
    return p_out
end

"""
Power of a polynomial.
"""
function ^(p::Polynomial, n::Int)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
    end
    return out
end

function ^(p::PolynomialModP, n::Int)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
    end
    return out
end

