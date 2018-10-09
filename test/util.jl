using Random


@testset "xlogx" begin
    @test iszero(xlogx(0))
    @test isnan(xlogx(NaN))

    @inferred xlogx(0)
    @inferred xlogx(0.)
    @inferred xlogx(0.0f0)
    
    for x = 1. : 1. : 10.
        @test xlogx(x) ≈ x * log(x)
    end
end


@testset "xlogy" begin
    @test iszero(xlogy(0, 1))
    @test iszero(xlogy(0, 0))
    @test iszero(xlogy(0, Inf))

    @test isnan(xlogy(NaN, 0))
    @test isnan(xlogy(NaN, 1))
    @test isnan(xlogy(0, NaN))

    @test_throws DomainError xlogy(0, -1)
    @test_throws DomainError xlogy(1, -1)
    @test_throws DomainError xlogy(NaN, -1)
    
    @inferred xlogy(0,0)
    @inferred xlogy(1,1)
    @inferred xlogy(1.,0.)

    @test xlogy(2, 3) ≈ 2.0 * log(3.0)

    for x = -10. : 1. : 10., y = 1. : 1. 10.
        @test xlogy(x,y) ≈ x * log(y)
    end
end


@testset "xexpx" begin
    @test iszero(xexpx(-Inf))
    @test isnan(xexpx(NaN))

    @inferred xexpx(0)
    @inferred xexpx(0.)
    @inferred xexpx(0.0f0)

    for x = -10. : 1. : 10.
        @test xexpx(x) ≈ x * exp(x)
    end
end


@testset "log1pexp" begin
    @test log1pexp(2.0)    ≈ log(1.0 + exp(2.0))
    @test log1pexp(-2.0)   ≈ log(1.0 + exp(-2.0))
    @test log1pexp(10000)  ≈ 10000.0
    @test log1pexp(-10000) ≈ 0.0
    @inferred log1pexp(0)
    @inferred log1pexp(0.)
    @inferred log1pexp(0.0f0)
end


@testset "log1mexp" begin
    @test log1mexp(-1.0)  ≈ log1p(- exp(-1.0))
    @test log1mexp(-10.0) ≈ log1p(- exp(-10.0))
end


@testset "Fermi Dirac binding probability" begin
    for φ = -10. : 1. : 10.
        p = fermi_dirac_prob(φ)
        
        @test p ≈ 1 / (1 + exp(-φ))

        @test fermi_dirac_logp(φ) ≈ log(p)
        @test fermi_dirac_1mp(φ) ≈ 1 - p
        @test fermi_dirac_l1mp(φ) ≈ log(1 - p)

        μ = randn()
        E = randn()
        φ = μ - E
        
        for f in (fermi_dirac_prob, 
                  fermi_dirac_logp, 
                  fermi_dirac_1mp, 
                  fermi_dirac_l1mp)
            
            @test f(μ, E) ≈ f(φ)        
        end
    end
end
