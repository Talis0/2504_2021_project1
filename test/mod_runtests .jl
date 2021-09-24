#############################################################################
#############################################################################
#
# A script that runs all unit tests in the project.
#                                                                               
#############################################################################
#############################################################################

include("../poly_factorization_project.jl")

####
# Execute unit tests for integers
###
include("mod_integers_test.jl")
mod_test_euclid_ints()
mod_test_ext_euclid_ints()

####
# Execute unit tests for polynomials
####
include("mod_polynomials_test.jl")
for i in [3,11,17,19,67,101]
    println("")
    println("testing p = $i")
    mod_sum_test_poly(i)
    mod_prod_test_poly(i)
    mod_prod_derivative_test_poly(i)
    mod_ext_euclid_test_poly(i)
    mod_division_test_poly(i)

end

####
# Execute unit tests for polynomial factorization
####
include("mod_factorization_test.jl")
mod_factor_test_poly()