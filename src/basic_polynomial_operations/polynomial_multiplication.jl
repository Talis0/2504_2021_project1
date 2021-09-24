#############################################################################
#############################################################################
#
# This file implements polynomial multiplication 
#                                                                               
#############################################################################
#############################################################################
using Primes

"""
Multiply two polynomials.
"""

#Multiplication in Z_p[x]

function *(p1::PolynomialModP, p2::PolynomialModP) #Fastest Method

    @assert p1.mod == p2.mod
    prime = p1.mod
    iszero(p1) && return zero(PolynomialModP,prime)
    iszero(p2) && return zero(PolynomialModP,prime)

    degree(p1) == 0 && return p2*leading(p1).coeff
    degree(p2) == 0 && return p1*leading(p2).coeff

    length(p1) == 0 && return p1.polynomial.terms[1]*p2
    length(p2) == 0 && return p2.polynomial.terms[1]*p1
    
    #Converts Polynomials into vector of coefficients
    p1 = coeffvector(p1.polynomial)
    p2 = coeffvector(p2.polynomial)

    n1 = length(p1)
    n2 = length(p2)
    output = zeros(n1+n2)
    outpoly = Term[]

    #Expands terms like with old multiplication method
    for i in 1:n1
        if p1[i] == 0
            i = i
        else 
            for j in 1:n2
                if p2[j] == 0
                    j = j
                else 
                    output[i+j] = mod(output[i+j] + mod(p1[i]*p2[j],prime),prime)
                end
            end
        end
    end
    output = Int.(output)

    #converting back to Polynomial
    for i in 1:n1+n2
        push!(outpoly, Term(output[i],(n1+n2)-i))
    end

    return PolynomialModP(outpoly,prime)

end 

function mult(p1::PolynomialModP, p2::PolynomialModP)::PolynomialModP #Old Method
    @assert p1.mod == p2.mod
    p = p1.mod
    p_out = PolynomialModP(p)
    for t in p1
        p_out = p_out + (t * p2)
    end
    return p_out
end  

function mult2(p1::PolynomialModP, p2::PolynomialModP) #Splitting method V2 (Demonstated in lectures)
    @assert p1.mod == p2.mod
    prime = p1.mod

    (p1 == 0 || p2 == 0) && return zero(prime)

    n1 = length(p1)
    n2 = length(p2)
    p1pairs = [(0,0) for _ in 1:n1]
    p2pairs = [(0,0) for _ in 1:n2]
    i = 1
    j = 1
    for t in p1.polynomial.terms
        p1pairs[i] = (t.coeff,t.degree)
        i = i+1
    end
    for t in p2.polynomial.terms
        p2pairs[j] = (t.coeff,t.degree)
        j = j+1
    end

    #Converts the polynomials to vetcors of tuples 

    return split(p1pairs,p2pairs,prime)
    
end

function mult3(p1::PolynomialModP, p2::PolynomialModP) #splitting method V1, same as mult 3 except with polynomials instead of vectors
    @assert p1.mod == p2.mod
    prime = p1.mod

    (p1 == 0 || p2 == 0) && return zero(prime)
    n1 = length(p1.polynomial.terms)
    n2 = length(p2.polynomial.terms)

    (n1 == 1 && n2 == 1) && return PolynomialModP(leading(p1)*leading(p2),prime)

    n1 == 1 && return p2*leading(p1)
    n2 == 1 && return p1*leading(p2)

    p1a = PolynomialModP(Polynomial(p1.polynomial.terms[1:Int(floor(n1/2))],true),prime,true)
    p1b = PolynomialModP(Polynomial(p1.polynomial.terms[Int(floor(n1/2+1)):n1],true),prime,true)
    p2a = PolynomialModP(Polynomial(p2.polynomial.terms[1:Int(floor(n2/2))],true),prime,true)
    p2b = PolynomialModP(Polynomial(p2.polynomial.terms[Int(floor(n2/2+1)):n2],true),prime,true)

    return (p1a*p2a+p1b*p2b+p1b*p2a+p2b*p1a)
    
end 

# Multiplication in Z[x]
height(p::Polynomial)::Int = maximum(abs.(coeff.(p.terms)))

nonzero(p::Polynomial)::Int = length(p.terms) #number of non-zero terms (the # operator)

required_M(p1::Polynomial, p2::Polynomial)::Int = 2 * height(p1) * height(p2) * minimum([nonzero(p1),nonzero(p2)]) + 1
#An upper bound on the product of primes used for CRT, M


function fullcrt(u::Vector{Int},m::Vector{Int}) #An implimentation of the full CRT algorthm
    n = length(u)
    v = Array{Int}(undef, n)
    v[1] = u[1]
    out = v[1]
    if n > 1
        v[2] = ((u[2]- v[1])*int_inverse_mod(m[1],m[2]))%m[2]
        modprod = m[1]
        a = v[1]
        out  = out+v[2]*modprod
        if n > 2
            for i in 2:n-1
                a = a + v[i]*modprod
                modprod = modprod*m[i]
                v[i+1] = ((u[i+1] - a)*int_inverse_mod(modprod,m[i+1]))%m[i+1]
                out  = out+v[i+1]*modprod
            end
        end
    end
    return symmod(out,modprod*last(m))
end

function crt(u::Vector{Int},m::Vector{Int}) #An implimentation of CRT for just two cases
    v = ((u[2]- u[1])*int_inverse_mod(m[1],m[2]))%m[2]
    out  = u[1]+v*m[1]
    return symmod(out,m[1]*m[2])
end

function split(a::Vector{Tuple{Int64, Int64}}, b::Vector{Tuple{Int64, Int64}}, p::Int)
    n1 = length(a)
    n2 = length(b)
    #used for mult2  and 3
    #Impliments the splitting multiplication demonstated in lecturs

    if (n1 == 1 && n2 == 1) #these are mostly edge cases 
    
        return PolynomialModP([Term(mod(a[1][1]*b[1][1],p),a[1][2]+b[1][2])],p)

    elseif (n1 == 0 || n2 == 0) 

         return zero(PolynomialModP,p)

    elseif (n1 == 1)
        

        coef = a[1][1] 
        deg = a[1][2]  
        i = 1
        for t in a
            t = (mod(t[1]*coef,p),t[2]+deg)
        end
        output = Term.(b)
        return PolynomialModP(output,p)

    elseif (n2 == 1)
 
        coef = b[1][1]  
        deg = b[1][2]  
        i = 1
        for t in b
            t = ( mod(t[1]*coef,p), t[2]+deg )
        end
        output = Term.(b)
        return PolynomialModP(output,p)

    else #splits vector in half and expands the terms uisng a recusive call
        a1 = a[1:Int(floor(n1/2))]
        a2 = a[Int(floor(n1/2+1)):n1]
        b1 = b[1:Int(floor(n2/2))]
        b2 = b[Int(floor(n2/2+1)):n2]

        return split(a1,b1,p) + split(a2,b2,p) + split(a1,b2,p) + split(a2,b1,p)
        
    end
end
function *(p1::Polynomial, p2::Polynomial) #CRT multiplication
    iszero(p1)||iszero(p2) && return zero(Polynomial)
    degree(p1) == 0 && return p2*leading(p1).coeff
    degree(p2) == 0 && return p1*leading(p2).coeff
    (degree(p1) + degree(p2) < 6) && (return p1*̄p2) #The old multimplication is more efficien for small degrees
    M = required_M(p1,p2)
    n = degree(p1)+degree(p2)
    prime = 1
    pprod = 3
    nextpolyterms = coeffvector(PolynomialModP(p1,3)*PolynomialModP(p2,3), n+1)

    while pprod < M

        prime = nextprime(pprod+1)
        nextpolyterms = mod.(nextpolyterms,pprod)

        p1mod = PolynomialModP(p1,prime)
        p2mod = PolynomialModP(p2,prime)
        p = p1mod*p2mod

        a = coeffvector(p,n+1)

        for j in 1:n+1
            t = crt([a[j],nextpolyterms[j]],[prime,pprod])
            nextpolyterms[j] = t
        end  
        pprod = prime*pprod
    end

    output = Array{Term}(undef,n+1)
    i = n+1

    for t in nextpolyterms
        output[i] = Term(t,i-1)
        i = i-1
    end

    return Polynomial(output)
end 

function *̄(p1::Polynomial, p2::Polynomial)::Polynomial #Old Multiplication
    p_out = Polynomial()
    for t in p1
        p_out = p_out + (t * p2)
    end
    return p_out
end

"""
Power of a polynomial.
"""

function pow(p::PolynomialModP, n::Int) #old ^ function, 
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
    end
    return out
end

function ^(p::Polynomial,n::Int) #The new optimised Power Function, which uses 
    n == 0 && return one(Polynomial) #the binary expansion of the power to reduce the number of operations needed
    n == 1 && return p
    iszero(p) && return zero(Polynomial)

    a = bitstring(n)
    twopowers = zeros(64)
    maxdegree = 1
    result = one(Polynomial)

    for i in 1:64
        (a[65-i] == '1') && (maxdegree = i)
        twopowers[i] = Int(a[i])-48        
    end
    #Makes a  vector that represents the binary expansion of the power

    (twopowers[64] == 1) && (result = result*p)

    for i in 2:maxdegree # loops through and multiplies together each of powers of two corresponding
                        # to the binary expansion
        p = p*p
        (twopowers[65-i] == 1) && (result = result*p)
    end

    return result 
end

function ^(p::PolynomialModP,n::Int) #Same as above implimentation accept with mod p 
    n == 0 && return one(PolynomialModP,p.mod)
    n == 1 && return p
    iszero(p) && return zero(PolynomialModP,p.mod)
    prime = p.mod
    a = bitstring(n)
    twopowers = zeros(64)
    maxdegree = 1
    for i in 1:64
 
        if (a[65-i] == '1') 
            maxdegree = i
            
        end
        twopowers[i] = Int(a[i])-48        
    end

    result = one(PolynomialModP,prime)

    for i in 1:maxdegree 
        (twopowers[65-i] == 1) && (result = result*p)
        p = p*p #everything gets reduced mod prime at this step

    end
    return result 
end



