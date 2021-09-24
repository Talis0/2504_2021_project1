#############################################################################
#############################################################################
#
# This file implements polynomial division 
#                                                                               
#############################################################################
#############################################################################

"""  Modular algorithm.
f divide by g

f = q*g + r

p is a prime
"""

function divide(num::Polynomial, den::Polynomial)
    function division_function(p::Int)
        f, g = mod(num,p), mod(den,p)
        degree(f) < degree(num) && return nothing 
        iszero(g) && throw(DivideError())
        q = Polynomial()
        prev_degree = degree(f)
        while degree(f) ≥ degree(g) 
            h = Polynomial( (leading(f) ÷ leading(g))(p) )  #syzergy 
            f = mod((f - h*g), p)
            q = mod((q + h), p)  
            prev_degree == degree(f) && break
            prev_degree = degree(f)
        end
        
        @assert iszero( mod((num  - (q*g + f)),p))
        return q, f
    end
    return division_function
end 

function divide(num::PolynomialModP, den::PolynomialModP)
    p = num.mod
    @assert num.mod == den.mod
    f, g = num, den
    if degree(g) > degree(f)
        output =  [zero(PolynomialModP,p) , f]
        return output
    end

    iszero(g) && throw(DivideError())
    iszero(f) && return zero(PolynomialModP,p), zero(PolynomialModP,p)
    

    q = PolynomialModP(f.mod)
    prev_degree = degree(f)
    while degree(f) ≥ degree(g) 
        h = PolynomialModP(Polynomial((leading(f) ÷ leading(g))(p)),p,true )  #syzergy 
        hg = - h*g
        f = plus(f,hg)
        q = q+h
        prev_degree == degree(f) && break
        prev_degree = degree(f)
    end

    @assert iszero(num - (q*g + f))
    
    return q,f 
end


"""
The quotient from polynomial division. Returns a function of an integer.
"""
÷(num::Polynomial, den::Polynomial)  = (p::Int) -> first(divide(num,den)(p)) 
÷(num::PolynomialModP, den::PolynomialModP) = first(divide(num,den))

"""
The remainder from polynomial division. Returns a function of an integer.
"""
rem(num::Polynomial, den::Polynomial)  = (p::Int) -> last(divide(num,den)(p))
rem(num::PolynomialModP, den::PolynomialModP)  = last(divide(num,den))