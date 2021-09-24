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
function mod_factor_test_poly(;N::Int = 25, seed::Int = 10, primes::Vector{Int} = [3])
    Random.seed!(seed)
    for prime in primes
        print("\ndoing prime = $prime \t")
        for _ in 1:N
            println(".")
            p = rand(PolynomialModP,prime)
            @show p
            factorization = factor(p)
            pr = expand_factorization(factorization)
            if p - pr != 0
                @assert p-pr == 0 
            end

        end
    end
    println("\nfactor_test_poly - PASSED")
end

mod_factor_test_poly()
