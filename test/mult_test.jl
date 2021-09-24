include("mod polynomials_test.jl")

function poly_benchmark()
    println("---- Multiplication Test ----")
    println("(Time for New method on to")
    println("")
    for mod in [5,41,211,1009,10007]
        println("Mod = $mod")
        println("")
        for i in [5,20,50,100,500]
            println("   Degree = $i:")
            print("     ")
            @time for _ in 1:50
                p1 = rand(PolynomialModP,mod, degree = i)
                p2 = rand(PolynomialModP,mod, degree = i)
                p1*p2
            end
            print("     ")
            @time for _ in 1:50
                p1 = rand(PolynomialModP,mod, degree = i)
                p2 = rand(PolynomialModP,mod, degree = i)
                p1*p2
            end
        end
        println("")
        @time mod_prod_test_poly(mod)
        println("")
    end
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

function big_poly_test(;n = 100::Int,coeff = 10000::Int)
    println("One Big Polynomial")
    p1 = rand(Polynomial, degree = n,max_coeff = coeff)
    p2 = rand(Polynomial, degree = n,max_coeff = coeff)
    @time a = p1*p2
    @time b = p1*̄p2
    @assert a == b
    
    println("")

end

function power_test(mod::Int, power::Int = 100)
    p1 = rand(PolynomialModP,mod)
    @assert p1^power == pow(p1,power)
    println("Power Test: PASSED")
end

#code profiler
