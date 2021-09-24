#############################################################################
#############################################################################
#
# This file contains units tests for polynomial operations
#                                                                               
#############################################################################
#############################################################################
"""
Test Polynomial Addition and Subtraction
"""

function mod_sum_test_poly(mod::Int ; num_sums::Int = 100, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:num_sums
        p1 = rand(PolynomialModP,mod)
        p2 = rand(PolynomialModP,mod)
        sum = p1+p2
        @assert p2 == sum - p1
        @assert p1 == sum - p2
    end
end
"""
Test product of polynomials.
"""
function mod_prod_test_poly(mod::Int ;N::Int = 100, N_prods::Int = 20, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialModP,mod)
        p2 = rand(PolynomialModP,mod)
        prod = p1*p2
        prodtest= leading(p1)*leading(p2)
        c = prodtest.coeff%mod
        prod_of_leading_terms = Term(c,prodtest.degree)
        @assert leading(prod) == prod_of_leading_terms
    end

    for _ in 1:N
        p_base = PolynomialModP(Term(1,0),mod)
        for _ in 1:N_prods
            p = rand(PolynomialModP,mod)
            prod = p_base*p
            prodtest= leading(p_base)*leading(p)
            c = prodtest.coeff%mod
            prodtest = Term(c,prodtest.degree)
            @assert leading(prod) == prodtest
            p_base = prod
        end
    end
    println("prod_test_poly - PASSED")
end

"""
Test derivative of polynomials (as well as product).
"""
function mod_prod_derivative_test_poly(mod::Int;N::Int = 10^2,  seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialModP,mod)
        p2 = rand(PolynomialModP,mod)
        p1d = derivative(p1)
        p2d = derivative(p2)
        @assert (p1d*p2) + (p1*p2d) == derivative(p1*p2)
    end
    println("prod_derivative_test_poly - PASSED")
end

"""
Test division of polynomials modulo p.
"""
function mod_division_test_poly(mod::Int; N::Int = 10^4, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialModP,mod)
        p2 = rand(PolynomialModP,mod)
        p_prod = p1*p2
        q, r = PolynomialModP(mod), PolynomialModP(mod)
        try
            q, r = divide(p_prod, p2)
            if (q, r) == (nothing,nothing)
                println("Unlucky prime: $p1 is reduced to $(p1 % prime) modulo $prime")
                continue
            end
        catch e
            if typeof(e) == DivideError
                @assert p2, prime == 0
            else
                throw(e)
            end
        end
        @assert iszero( q*p2+r - p_prod )
    end
    println("division_test_poly - PASSED")
end

"""
Test the extended euclid algorithm for polynomials modulo p.
"""
function mod_ext_euclid_test_poly(mod::Int; N::Int = 100, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialModP,mod)
        p2 = rand(PolynomialModP,mod)
        g, s, t = extended_euclid_alg(p1, p2)
        @assert s*p1 + t*p2 - g == 0
    end
    println("ext_euclid_test_poly - PASSED")
end