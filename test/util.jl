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


@testset "Fermi Dirac binding probability" begin
    Random.seed!(119808)
    
    for _ = 1 : 10
        φ = randn()
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
