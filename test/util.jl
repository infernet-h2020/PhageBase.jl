using Random


@testset "xlogx, xexpx" begin
    @test iszero(xlogx(0))

    @test xlogx(0) isa Float64
    @test xlogx(Float16(0)) isa Float16

    @test iszero(xexpx(-Inf))

    @test xexpx(0) isa Float64
    @test xexpx(Float16(0)) isa Float16
    
    for testid = 1:10
        Random.seed!(testid + 160)
        x = 10rand()
        iszero(x) && continue   # should never happen
        @test xlogx(x) ≈ x * log(x)
        @test xexpx(x) ≈ x * exp(x)
    end
end


@testset "log1pexp" begin

    @test log1pexp(1) isa Float64
    @test log1pexp(Float16(1)) isa Float16
    @test log1pexp(Int(1)) isa Float64

    @test log1pexp(1) ≈ 1.31326168751822283404899549496785564191528008567034837471
    @test log1pexp(-1) ≈ 0.313261687518222834048995494967855641915280085670348374719
    
    # test all branches
    @test log1pexp(-40) ≈ 4.24835425529158898630497784363158218187775067982297454220e-18
    @test log1pexp(0) ≈ 0.69314718055994530941723212145817656807550013436025525412068000
    @test log1pexp(20) ≈ 20.000000002061153620314380703238982798877915235602669875387647193
    @test log1pexp(40) ≈ 40.0000000000000000042483542552915889863049778436315821818777506798
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
