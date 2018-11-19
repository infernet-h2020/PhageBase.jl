@testset "prior" begin
    for testrep=1:5
        Random.seed!(2176+testrep)
        A=rand(2:6); L=rand(2:6);
        flen = fieldslen(A,L)

        η = rand()

        prior = GaussPrior{A,L}(η)
        @test all(x -> x == η, prior.η)
        @test length(prior.η) == flen
        
        randn!(prior.η)
        randn!(prior.ξ)

        fields = Fields{A,L}(randn(flen + 1))
        @test log_prior(fields, prior) ≈ -0.5sum(prior.η .* (fields.x[1:flen] .- prior.ξ) .* (fields.x[1:flen] .- prior.ξ))

        f(x) = log_prior(Fields{A,L}(x), prior)
        Gad = ForwardDiff.gradient(f, fields.x)
        G = rand(flen + 1)
        G0 = copy(G)
        log_prior_grad!(G, fields, prior)

        @test G .- G0 ≈ Gad
        @test iszero(Gad[end])

        Had = ForwardDiff.hessian(f, fields.x)
        H = rand(flen + 1, flen + 1)
        H0 = copy(H)
        log_prior_hess!(H, prior)

        @test (H .- H0) ≈ Had
        @test iszero(Had[end,:]) && iszero(Had[:,end])
    end
end
