#############################################################################
#############################################################################
#
# This file defines the polynomial type with several operations 
#                                                                               
#############################################################################
#############################################################################

####################################
# Polynomial and PolynomialModP Type #
#################################### 

    """
    A Polynomial type - designed to be for polynomials with integer coefficients.
    """
    struct Polynomial
        terms::Array{Term}   
        Polynomial() = new(Term[])
        #Inner constructor
        Polynomial(h::Array{Term}, clean::Bool) = clean && new(h) 
    end

    """
    Construct a polynomial with a single term.
    """
    function Polynomial(t::Term)
        terms = Term[]
        t.coeff != 0 && push!(terms, t)
        return Polynomial(terms,true)
    end

    """
    Construct a polynomial with a vector of terms.
    """
    function Polynomial(tv::Vector{Term})
        n = length(tv)
       sort!(tv, by = x -> -x.degree)
        if n == 1
            return Polynomial(tv, true)
        else
            terms = Term[]
            i = 1
            while i < n+1
                if tv[i].coeff != 0 
                    if i == n
                        push!(terms,tv[i])

                    elseif tv[i].degree == tv[i+1].degree 
                        tv[i+1] = tv[i]+tv[i+1]
                    else
                        push!(terms,tv[i])
                    end
                end
                i = i + 1
            end
        end
        
        return Polynomial(terms, true) 
        #Once made pretty (no repreated factors and ordered), the vector is sent to the inner constructor to 
        #made into a polynomial
    end 

    """
    Polymonials mod P
    """
    struct PolynomialModP 
        polynomial::Polynomial
        mod::Int
        PolynomialModP(p) = new(Polynomial(Term[]),p)
        PolynomialModP(poly::Polynomial, mod::Int, modded::Bool) = modded && new(poly,mod)
    end

    function PolynomialModP(t::Term,p)
        terms = Term[]
        modcoeff = mod(t.coeff,p)
        modcoeff != 0 && push!(terms,Term(modcoeff,t.degree))
        return PolynomialModP(Polynomial(terms),p,true)
    end

    function PolynomialModP(tv::Vector{Term},p::Int)
        n = length(tv)
        sort!(tv, by = x -> -x.degree)
        tv = mod.(tv,p)
        if n == 1
            terms = tv
        else
            terms = Term[]
            i = 1
            while i < n+1
                if tv[i].coeff != 0 
                    tv[i].coeff 
                    if i == n
                        push!(terms,tv[i])
                    elseif tv[i].degree == tv[i+1].degree 
                        tv[i+1] = tv[i]+tv[i+1]
                        tv[i+1].coeff = tv[i+1].coeff%p
                    else
                        push!(terms,tv[i])
                    end
                end
                i = i + 1
            end
        end

        if length(terms)== 0 || terms[1].coeff == 0 
            return PolynomialModP(p)
        end
        
        terms = Polynomial(terms,true)
        return PolynomialModP(terms, p, true) 
    end 

    PolynomialModP(poly::Polynomial, p::Int) = PolynomialModP(mod(poly,p), p, true) 

####################################
# Polynomial Construction #
#################################### 
    """
    Construct a polynomial of the form x^p-x.
    """
        cyclotonic_polynomial(p::Int) = Polynomial([Term(1,p), Term(-1,0)])

    """
    Construct a polynomial of the form x-n.
    """
        linear_monic_polynomial(n::Int) = Polynomial([Term(1,1), Term(-n,0)])

    """
    Construct a polynomial of the form x.
    """
    x_poly() = Polynomial(Term(1,1))
    x_poly(p::Int) = PolynomialModP(Term(1,1),p)

    """
    Creates the zero polynomial.
    """
    zero(::Type{Polynomial})::Polynomial = Polynomial()
    zero(::Type{PolynomialModP},n::Int)::PolynomialModP = PolynomialModP(n)
    """
    Creates the unit polynomial.
    """
    one(::Type{Polynomial})::Polynomial = Polynomial(Term(1,0))
    one(p::Polynomial) = one(typeof(p))

    one(::Type{PolynomialModP},p)::PolynomialModP = PolynomialModP(one(Term),p)
    one(p::PolynomialModP) = one(typeof(p),p.mod)

    """
    Generates a random polynomial.
    """
    function rand(::Type{Polynomial} ; 
                    degree::Int = -1, 
                    terms::Int = -1, 
                    max_coeff::Int = 100, 
                    mean_degree::Float64 = 5.0,
                    prob_term::Float64  = 0.7,
                    monic = false,
                    condition = (p)->true)
            
        while true 
            _degree = degree == -1 ? rand(Poisson(mean_degree)) : degree
            _terms = terms == -1 ? rand(Binomial(_degree,prob_term)) : terms
            degrees = vcat(sort(sample(0:_degree-1,_terms,replace = false)),_degree)
            coeffs = rand(1:max_coeff,_terms+1)
            monic && (coeffs[end] = 1)
            p = Polynomial( [Term(coeffs[i],degrees[i]) for i in 1:length(degrees)] )
            condition(p) && return p
        end
    end

    function rand(::Type{PolynomialModP}, mod::Int ; 
        degree::Int = -1, 
        terms::Int = -1, 
        mean_degree::Float64 = 5.0,
        prob_term::Float64  = 0.7,
        monic = false,
        condition = (p)->true)

        while true 
            max_coeff = mod
            _degree = degree == -1 ? rand(Poisson(mean_degree)) : degree
            _terms = terms == -1 ? rand(Binomial(_degree,prob_term)) : terms
            degrees = vcat(sort(sample(0:_degree-1,_terms,replace = false)),_degree)
            coeffs = rand(1:max_coeff-1,_terms+1)
            monic && (coeffs[end] = 1)
            p = [Term(coeffs[i],degrees[i]) for i in 1:length(degrees)]
            poly = PolynomialModP(p,mod)
            condition(poly) && return poly
        end
    end

###########
# Display #
###########

    """
    Show a polynomial.
    """
    function show(io::IO, p::Polynomial) 
        p = deepcopy(p)
        
        if iszero(p)
            print(io,"0")
        else
            n = length(p.terms)
            for (i,t) in enumerate(p.terms)
            
                if t.coeff < 0
                    print(io,i != 1 ? " - " : "-")
                    print(io,-t)
                else 
                    print(io,i != 1 ? " + " : "")
                    print(io,t)
                end
            end
        end
    end

    show(io::IO,poly::PolynomialModP) = print(io, poly.polynomial, " (mod ", poly.mod, ")")


##############################################
# Iteration over the terms of the polynomial #
##############################################

    """
    Allows to do iteration over the terms of the polynomial. The iteration is in an arbitrary order.
    """
    iterate(p::Polynomial, state=1) = iterate(p.terms, state)
    iterate(p::PolynomialModP, state=1) = iterate(p.polynomial.terms, state)

    function coeffvector(p::Polynomial, size::Int = -1)::Vector{Int}
        size == -1 && (size = degree(p)+1)
        a = zeros(size)
        if length(p.terms) != size
            for t in p.terms
                a[size - t.degree] = t.coeff
            end
        else
            a = coeffs(p)
        end
        return a
    end

    coeffvector(p::PolynomialModP,size::Int = -1) = coeffvector(p.polynomial, size)
 
##############################
# Queries about a polynomial #
##############################

    """
    The number of terms of the polynomial.
    """
    length(p::Polynomial) = length(p.terms)
    length(p::PolynomialModP) = length(p.polynomial)

    """
    The leading term of the polynomial.
    """
    leading(p::Polynomial)::Term = isempty(p.terms) ? zero(Term) : first(p.terms) 
    leading(p::PolynomialModP)::Term = leading(p.polynomial) 

    """
    Returns the coefficients of the polynomial.
    """
    coeffs(p::Polynomial)::Vector{Int} = [t.coeff for t in p]
    coeffs(p::PolynomialModP)::Vector{Int} = coeffs(p.polynomial)


    """
    The degree of the polynomial.
    """
    degree(p::Polynomial)::Int = leading(p).degree 
    degree(p::PolynomialModP)::Int = leading(p.polynomial).degree 

    """
    The content of the polynomial is the GCD of its coefficients.
    """
    content(p::Polynomial)::Int = euclid_alg(coeffs(p))
    content(p::PolynomialModP)::Int = content(p.polynomial) 

    """
    Evaluate the polynomial at a point `x`.
    """
    evaluate(f::Polynomial, x::T) where T <: Number = sum(evaluate(t,x) for t in f)
    evaluate(f::PolynomialModP, x::T) where T <: Number = evaluate(f.polynomial,x)%f.mod

################################
# Pushing and popping of terms #
################################

    """
    Push a new term into the polynomial.
    """
    #Note that ideally this would throw and error if pushing another term of degree that is already in the polynomial
    function push!(p::Polynomial, t::Term) 
        iszero(t) && return #don't push a zero
        Polynomial(push!(p.terms,t))
    end

    function push!(p::PolynomialModP, t::Term)
        iszero(t) && return #don't push a zero
        PolynomialModP(push!(p.polynomial,t),p.mod)
    end

    """
    Pop the leading term out of the polynomial.
    """
    pop!(p::Polynomial)::Term = pop!(p.terms)
    pop!(p::PolynomialModP)::Term = pop!(p.polynomial.terms)

    prepop!(p::Polynomial)::Term = popat!(p.terms,1) #I needed a way to "pop" the first term.
    prepop!(p::PolynomialModP)::Term = popat!(p.polynomial.terms,1)
    """
    Check if the polynomial is zero.
    """
    iszero(p::Polynomial)::Bool = isempty(p.terms)
    iszero(p::PolynomialModP)::Bool = iszero(p.polynomial)

#################################################################
# Transformation of the polynomial to create another polynomial #
#################################################################

    """
    The negative of a polynomial.
    """
    -(p::Polynomial) = Polynomial(map((pt)->-pt, p.terms))
    -(p::PolynomialModP) = PolynomialModP(-p.polynomial,p.mod)

    """
    Create a new polynomial which is the derivative of the polynomial.
    """
    function derivative(p::Polynomial)::Polynomial 
        der_p = Polynomial()
        for term in p
            der_term = derivative(term)
            !iszero(der_term) && push!(der_p,der_term)
        end
        return der_p
    end

    derivative(p::PolynomialModP)::PolynomialModP = PolynomialModP(derivative(p.polynomial),p.mod)

    """
    The prim part (multiply a polynomial by the inverse of its content).
    """
    prim_part(p::Polynomial) = p ÷ content(p)
    prim_part(p::PolynomialModP) = p ÷ content(p)


    """
    A square free polynomial.
    """
    square_free(p::Polynomial, prime::Int)::Polynomial = (p ÷ gcd(p,derivative(p),prime))(prime)
    square_free(p::PolynomialModP)::PolynomialModP = p ÷ gcd(p,derivative(p))

#################################
# Queries about two polynomials #
#################################

    """
    Check if two polynomials are the same
    """
    ==(p1::Polynomial, p2::Polynomial)::Bool = p1.terms == p2.terms
    ==(p1::PolynomialModP, p2::PolynomialModP)::Bool = p1.polynomial.terms == p2.polynomial.terms


    """
    Check if a polynomial is equal to 0.
    """
    #Note that in principle there is a problem here. E.g The polynomial 3 will return true to equalling the integer 2.
    ==(p::Polynomial, n::T) where T <: Real = iszero(p) == iszero(n)
    ==(p::PolynomialModP, n::T) where T <: Real = ==(p.polynomial,n)


##################################################################
# Operations with two objects where at least one is a polynomial #
##################################################################
    """
    Subtraction of two polynomials.
    """
    -(p1::Polynomial, p2::Polynomial)::Polynomial = p1 + (-p2)
    -(p1::PolynomialModP, p2::PolynomialModP)::PolynomialModP = p1 + (-p2)

    -(p::PolynomialModP, n::Int) = p + PolynomialModP(Term(-n,0),p.mod) 
    #Subtract an integer from a Polynomial
    -(p::Polynomial, n::Int)::Polynomial = p + Polynomial(Term(-n,0)) 


    """
    Multiplication of polynomial and term.
    """
    *(t::Term,p1::Polynomial)::Polynomial = iszero(t) ? Polynomial() : Polynomial(map((pt)->t*pt, p1.terms))
    *(p1::Polynomial, t::Term)::Polynomial = t*p1

    *(t::Term,p1::PolynomialModP) = iszero(t) ? PolynomialModP(p1.mod) : PolynomialModP(map((pt)->mod(t*pt,p1.mod), p1.polynomial.terms),p1.mod)
    *(p1::PolynomialModP, t::Term) = t*p1

    """
    Multiplication of polynomial and an integer.
    """
    *(n::Int,p::Polynomial)::Polynomial = p*Term(n,0)
    *(p::Polynomial,n::Int)::Polynomial = n*p

    *(n::Int,p::PolynomialModP) = p*Term(n,0)
    *(p::PolynomialModP,n::Int) = n*p

    """
    Integer division of a polynomial by an integer.

    Warning this may not make sense if n does not divide all the coefficients of p.
    """
    ÷(p::Polynomial,n::Int) = (prime)->Polynomial(map((pt)->((pt ÷ n)(prime)), p.terms))
    ÷(p::PolynomialModP, n::Int) = PolynomialModP(÷(p.polynomial,n)(p.mod),p.mod)

    """
    Take the s of a polynomial with an integer.
    """
    function mod(f::Polynomial, p::Int)::Polynomial
        p_out = Polynomial()
        for t in f
            push!(p_out, mod(t, p)) #if coeff reduced to zero, push! will handle it
        end
        return p_out
    end

    """
    Power of a polynomial mod prime.
    """
    function pow_mod(p::Polynomial, n::Int, prime::Int)
        n < 0 && error("No negative power")
        out = one(p)
        for _ in 1:n
            out *= p
            out = mod(out, prime)
        end
        return out
    end 

" ---- Polynomial Run ----"