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
function mod_factor_test_poly(;N::Int = 24, seed::Int = 0, primes::Vector{Int} = [3,11,101])
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
