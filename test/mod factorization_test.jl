#############################################################################
#############################################################################
#
# This file contains units tests for polynomial factorization
#                                                                               
#############################################################################
#############################################################################


"""
Test factorization of polynomials.
"""
function mod_factor_test_poly(;N::Int = 25, seed::Int = 2, primes::Vector{Int} = [5,7,11,17])
    Random.seed!(seed)
    for prime in primes
        print("\ndoing prime = $prime \t")
        for _ in 1:N
            print(".")
            p = rand(PolynomialModP,prime)
            factorization = factor(p)
            pr = expand_factorization(factorization)
            @assert p-pr == 0 
        end
    end

    println("\nfactor_test_poly - PASSED")
end

mod_factor_test_poly()
