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
        @assert iszero((p2+p1) - sum)
    end
    println("Addition Test - PASSED")
end
"""
Test product of polynomials.
"""
function mod_prod_test_poly(m::Int ;N::Int = 100, N_prods::Int = 20, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialModP,m)
        p2 = rand(PolynomialModP,m)
        prod = p1*p2
        prodtest= mod(leading(p1)*leading(p2),m)
        if !iszero(prodtest - leading(prod))
            println("There is a multiplication problem for $p1 & $p2")
            @assert iszero(prodtest - leading(prod))
        end
    end

    for _ in 1:N
        p_base = PolynomialModP(Term(1,0),m)
        for _ in 1:N_prods
            p = rand(PolynomialModP,m)
            prod = p_base*p
            test = leading(p_base)*leading(p)
            test.coeff = (test.coeff)%m
            if !iszero(test - leading(prod))
                println("There is a multiplication problem for $p1 & $p2")
                @assert iszero(test - leading(prod))
            end
            p_base = prod
        end
    end
    println("Multiplication Test - PASSED")
end

function new_prod_test(p1,p2)
    println("")
    a = p1*p2 - p1*̄p1
    if a != 0
        println("The multiplication does not work for:")
        println(p1*p2)
        println(p1*̄p2) 
        return 
    end
    return "Multiplication Works"
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
        @assert iszero((p1d*p2) + (p1*p2d)-derivative(p1*p2))
    end
    println("Derivative Test - PASSED")
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
    println("Division Test - PASSED")
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
    println("Extended Euclid Algorithm Test - PASSED")
end