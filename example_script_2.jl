include("poly_factorization_project.jl")
println("Welcome to the Great Polynomial Factoriser")
println("")
println("")
println("----- Example Polynomials -----")
x = x_poly()
p1 = 2x^3 + 4x^2 - 3x
p2 = -6x^4 - 4x^2 - 3x + 3
println(p1)
println(p2)
println("")
println("mod 7")
x = x_poly(7)
p1 = 2x^5 + 5x - 3
p2 = -6x^4 - 4x^2 - 3x + 3
println(p1)
println(p2)
println("")
println("mod 17")
x = x_poly(17)
p1 = 16x^3 - 3
p2 = -6x^8 - 4x^2 - 3x + 1
println(p1)
println(p2)
println("")
println("mod 101")
x = x_poly(101)
a = 99x^8+ 50x^2 - 3
b = 3x^7+ 5x^7 + 10x^2 + 60x 
@show a
@show b

println("")

println("----- Basic Operations -----")
@show a+b
@show a-b
@show a*b
@show a÷b
@show a^3
println("")
println("d/dx(", "$a) = ", derivative(a))
println("")
println("-----Factorisation Example------")
p = (7x^3 + 2x^2 + 8x + 1)*(x^2+x+1)

println("The Factorisation of: ", p," is:")
println("")
factorization = factor(p,true)
println("Here it is:")
pretty_factor(factorization)
pr = expand_factorization(factorization)
println("")
println("Reconstructing: ", pr)

println(" ")
println("---- Extended Euclid Algorithm on x^2 + 50 and 10x^3+10 (mod 101) ----")
egcd = extended_euclid_alg(x^2 + 50,10x^3+10)
pretty_print_egcd((((x^2 + 50)*(10x^3+10)).polynomial,b.polynomial), egcd)

println("")
println("----- Bechmarking -----")
println("Time taken to raise $a to the 1000th power:")
@time l = a^1000
println("")
multbenchZ()
println("")
multbenchZp()
println("")
powbench()
println("")

println(":)")