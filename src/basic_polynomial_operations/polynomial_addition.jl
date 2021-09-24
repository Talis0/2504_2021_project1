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

function plus(p1::Polynomial, p2::Polynomial)::Polynomial
    iszero(p1) && return p2
    iszero(p2) && return p1
    (iszero(p1) && iszero(p1)) && return p2

    degree(p1) > degree(p2) ? n = degree(p1)+1 : n = degree(p2)+1

    p1 = coeffvector(p1,n)
    p2 = coeffvector(p2,n)

    output = [Term(0,0) for _ in 1:n]
    for i in 1:n
        output[i] = Term(p1[i]+p2[i],n-i)
    end

    return Polynomial(output,true)

end

function plus(p1::PolynomialModP, p2::PolynomialModP)::PolynomialModP #old, slow addition method
    @assert p1.mod == p2.mod
    p = p1.mod
    p1, p2 = deepcopy(p1), deepcopy(p2)
    p3 = Term[]
    while !iszero(p1) && !iszero(p2)
        t1, t2 = leading(p1), leading(p2) 
        if t1.degree == t2.degree
            push!(p3, mod((prepop!(p1)+prepop!(p2)),p))
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
    
    return PolynomialModP(p3,p)
end

function +(p1::PolynomialModP, p2::PolynomialModP)::PolynomialModP #New, fast addition method
    @assert p1.mod == p2.mod
    p = p1.mod
    iszero(p1) && return p2
    iszero(p2) && return p1
    (iszero(p1) && iszero(p1)) && return p2

    degree(p1) > degree(p2) ? n = degree(p1)+1 : n = degree(p2)+1

    p1 = coeffvector(p1.polynomial,n)
    p2 = coeffvector(p2.polynomial,n)

    output = [Term(0,0) for _ in 1:n]
    i = 1
    for j in 1:n
        c = mod(p1[j]+p2[j],p)
        if c == 0
            popat!(output,i)
        else
            output[i] = Term(c,n-j)
            i = i+1
        end
        
    end

    poly = Polynomial(output,true)
    return PolynomialModP(poly,p,true)

end

"""
Add a polynomial and a term.
"""
+(p::Polynomial, t::Term) = p + Polynomial(t)
+(t::Term, p::Polynomial) = p + t

+(p::PolynomialModP, t::Term) = p + PolynomialModP(t,p.mod)
+(t::Term, p::PolynomialModP) = p.polynomial + t

"""
Add a polynomial and an integer.
"""
+(p::Polynomial, n::Int) = p + Term(n,0)
+(n::Int, p::Polynomial) = p + Term(n,0)

+(p::PolynomialModP, n::Int) = p + Term(n,0)
+(n::Int, p::PolynomialModP) = p + Term(n,0)


function addtest()
    a = rand(PolynomialModP,101, degree = 10000)
    b = rand(PolynomialModP,101, degree = 10000)
    println("Time taken for:")
    println("New Method:")
    @time l = b+a
    println("Old Method:")
    @time m = plus(b,a)
    @assert iszero(l-m)
end