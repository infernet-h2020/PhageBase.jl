using Random


@testset "xlogx, xexpx" begin

    @test iszero(xexpx(-Inf))

    @test xexpx(0) isa Float64
    @test xexpx(Float16(0)) isa Float16

    @inferred xexpx(-Inf)
    @inferred xexpx(0)
    @inferred xexpx(0.)
    @inferred xexpx(NaN)
    
    for testid = 1:10
        Random.seed!(testid + 160)
        x = 10rand()
        
        iszero(x) && continue   # should never happen

        @test xexpx(x) ≈ x * exp(x)
    end
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
