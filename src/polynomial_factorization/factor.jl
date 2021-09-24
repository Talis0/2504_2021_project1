#############################################################################
#############################################################################
#
# This file implements factorization 
#                                                                               
#############################################################################
#############################################################################

"""
Factors a polynomial over the field Z_p.

Returns a vector of tuples of (irreducible polynomials (mod p), multiplicity) such that their product of the list (mod p) is f. Irreducibles are fixed points on the function factor.
"""

function factor(f::PolynomialModP, show::Bool = false)::Vector{Tuple{PolynomialModP,Int}}
    
    #Cantor Zassenhaus factorization

    degree(f) ≤ 1 && return [(f,1)]

    # make f primitive
    ff = prim_part(f)     

     # make f square-free
    squares_poly = gcd(f, derivative(ff)) 
  
    ff = ff ÷ squares_poly
    show && println(" The No_Square Polynomial Based of f:")
    show && println(ff)
    show && println("")

    # make f monic
    old_coeff = leading(ff).coeff
    ff = ff ÷ old_coeff  
    
    show && println(" The monic Polynomial Based of f:")
    show && println(ff)
    show && println("")
  
    dds = dd_factor(ff)

    ret_val = Tuple{PolynomialModP,Int}[]

    for (k,dd) in enumerate(dds)
        sp = dd_split(dd, k)
        sp = map((p)->p ÷ leading(p).coeff,sp) #makes the polynomials inside the list sp, monic
        for mp in sp
            push!(ret_val, (mp, multiplicity(f,mp)) )
        end
    end

    #Append the leading coefficient as well
    push!(ret_val, (leading(f).coeff* one(PolynomialModP,f.mod), 1) )

    degreetest = 0
    for i in ret_val
        degreetest = degreetest + degree(i[1]) * i[2]
    end

    if degreetest != degree(f) #when there is a factor with a multiplicity of p (the modulus)
                         #it does not show up in the factorisation since its derivative is 0
                         #This tries find these lost factors.
        l = pop!(ret_val)
        temp_poly = 2*l[1]
        push!(ret_val,(temp_poly,1))

        d = degree(f) - degreetest #this is the degree of this missing polynomial
        x = x_poly(f.mod)
        p = f.mod 

        miss_poly =  f ÷ expand_factorization(ret_val) #Figuring out what the missing polynomial is so we can factorise it
  
        for i in 1:d÷p
            cand = gcd(x^(p^(i))-x,miss_poly) #looks for candidate factor on irriducable polys
            if degree(cand) == i  #if it is infact irriducable, we can use it straight away and add it to the factorisation
                
                push!(ret_val,(cand,p))

            elseif degree(cand) > 0 #else, it is the productt of multiple irriducable polys, ad thus we must factorise again

                a = factor(cand)
                for j in a

                    j[2] = j[2]*p
                end
                ret_val = vcat(ret_val,a)
            end
            
        end 
    end 

    return ret_val
end 


"""
Expand a factorization.
"""

function expand_factorization(factorization::Vector{Tuple{PolynomialModP,Int}})::PolynomialModP 
    length(factorization) == 1 && return first(factorization[1])^last(factorization[1])
    return *([first(tt)^last(tt) for tt in factorization]...)
end

"""
Compute the number of times g divides f
"""

function multiplicity(f::PolynomialModP, g::PolynomialModP)::Int
    (iszero(f) || iszero(g)) && return 0
    degree(gcd(f, g)) == 0 && return 0
    return 1 + multiplicity(f ÷ g, g)
end


"""
Distinct degree factorization.

Given a square free polynomial `f` returns a list, `g` such that `g[k]` is a product of irreducible polynomials of degree `k` for `k` in 1,...,degree(f) ÷ 2, such that the product of the list (mod `prime`) is equal to `f` (mod `prime`).
"""
function dd_factor(f::PolynomialModP)::Array{PolynomialModP}
    x = x_poly(f.mod)
    w = deepcopy(x)
    g = Array{PolynomialModP}(undef,degree(f.polynomial)) #Array of polynomials indexed by degree
    #Looping over degrees
    for k in 1:degree(f)
        w = rem(w^(f.mod), f)
        g[k] = gcd(w - x, f) 

        f = f ÷ g[k]
    end
    
    #edge case for final factor
    f != one(PolynomialModP,f.mod) && push!(g,f)

    return g
end

"""
Distinct degree split.

Returns a list of irreducible polynomials of degree `d` so that the product of that list (mod prime) is the polynomial `f`.
"""

function dd_split(f::PolynomialModP, d::Int)::Vector{PolynomialModP}

    degree(f) == d && return [f]
    degree(f) == 0 && return []

    w = rand(PolynomialModP, f.mod, degree = d, monic = true)

    n_power = ((f.mod)^d-1) ÷ 2

    g = gcd(w^n_power - one(PolynomialModP,f.mod), f)

    ḡ = f ÷ g # g\bar + [TAB]
   
    return vcat(dd_split(g, d), dd_split(ḡ, d) )
end