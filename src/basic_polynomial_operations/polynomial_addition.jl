#############################################################################
#############################################################################
#
# This file implements polynomial addition 
#                                                                               
#############################################################################
#############################################################################

"""
Add two polynomials.
"""
function +(p1::Polynomial, p2::Polynomial)::Polynomial
    p1, p2 = deepcopy(p1), deepcopy(p2)
    p3 = Term[]
    while !iszero(p1) && !iszero(p2)
        t1, t2 = leading(p1), leading(p2) 
        if t1.degree == t2.degree
            push!(p3, (prepop!(p1)+prepop!(p2)))
        elseif t1.degree < t2.degree
            push!(p3,pop!(p2))
        else
            push!(p3,pop!(p1))
        end
    end
    while !iszero(p1)
        push!(p3,pop!(p1))
    end
    while !iszero(p2)
        push!(p3,pop!(p2))
    end
    
    return Polynomial(p3)
end

function +(p1::PolynomialModP, p2::PolynomialModP)::PolynomialModP
    @assert p1.mod == p2.mod
    return PolynomialModP(p1.polynomial + p2.polynomial,p1.mod)
end

"""
Add a polynomial and a term.
"""
+(p::Polynomial, t::Term) = p + Polynomial(t)
+(t::Term, p::Polynomial) = p + t

+(p::PolynomialModP, t::Term) = PolynomialModP(p.polynomial + Polynomial(t),p.mod)
+(t::Term, p::PolynomialModP) = p.polynomial + t

"""
Add a polynomial and an integer.
"""
+(p::Polynomial, n::Int) = p + Term(n,0)
+(n::Int, p::Polynomial) = p + Term(n,0)

+(p::PolynomialModP, n::Int) = p + Term(n,0)
+(n::Int, p::PolynomialModP) = p + Term(n,0)
