include("mod_polynomials_test.jl")

function poly_mult_benchmark()
    println("---- Multiplication Test ----")
    println("(Time for New method on top")
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

function big_poly_Benchmark(;n = 100::Int,coeff = 10000::Int)
    println("Time Test for One Big Polynomial")
    p1 = rand(Polynomial, degree = n,max_coeff = coeff)
    p2 = rand(Polynomial, degree = n,max_coeff = coeff)
    @time a = p1*p2
    @time b = p1*̄p2
    @assert a == b
    println("")
end


function addbench()
    println("---- Time Compairson for the Two Addition Algorithms ----")
    a = rand(PolynomialModP,101, degree = 10000)
    b = rand(PolynomialModP,101, degree = 10000)

    println("")
    println("+ function:")
    @time for _ in 1:10000
        p1 = rand(PolynomialModP,101)
        p2 = rand(PolynomialModP,101)
        p1+p2
    end

    println("")
    println("Plus() function:")
    @time for _ in 1:10000
        p1 = rand(PolynomialModP,101)
        p2 = rand(PolynomialModP,101)
        plus(p1,p2)
    end
    
end

function multbenchZ()
    println("---- Time Compairson for the Two Multiplication Algorithms (Z[x]) ----")
    println("")
    println("CRT: *")
    @time for _ in 1:100
        p1 = rand(Polynomial,degree = 100)
        p2 = rand(Polynomial,degree = 100)
        p1*p2
    end
    println("")
    println("Expansion: *̄()")

    @time for _ in 1:100
        p1 = rand(Polynomial, degree = 100)
        p2 = rand(Polynomial,degree = 100)
        p1 *̄ p2
    end
    
end

function multbenchZp()
    println("---- Time Compairson for the Four Multiplication Algorithms ----")
    println("")
    println("Vector Distribution: *")
    @time for _ in 1:100
        p1 = rand(PolynomialModP,101,degree = 100)
        p2 = rand(PolynomialModP,101,degree = 100)
        p1*p2
    end
    println("")
    println("Polynomial Distribution: mult()")

    @time for _ in 1:100
        p1 = rand(PolynomialModP,101,degree = 100)
        p2 = rand(PolynomialModP,101,degree = 100)
        mult(p1,p2)
    end
    println("")
    println("Vector Splitting: mult2()")
    @time for _ in 1:100
        p1 = rand(PolynomialModP,101,degree = 100)
        p2 = rand(PolynomialModP,101,degree = 100)
        mult2(p1,p2)
    end
    println("")
    println("Polynomial Splitting: mult3()")
    @time for _ in 1:100
        p1 = rand(PolynomialModP,101,degree = 100)
        p2 = rand(PolynomialModP,101,degree = 100)
        mult3(p1,p2)
    end
end

function powbench()
    println("")
    println("---- Time Compairson for the two Power Algorithms Z_p[x]----")
    println("Binary Expansion: ^()")
    @time for _ in 1:10
        p1 = rand(PolynomialModP,101,degree = 10)
        p1^200
    end
    println("")
    println("Repeated Multiplication: pow()")
    @time for _ in 1:10
        p1 = rand(PolynomialModP,101,degree = 10)
        pow(p1,200)
    end
end