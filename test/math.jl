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

    for x = -10:10, y = 1:10
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


@testset "xexpy" begin
    @test iszero(xexpy(Inf, -Inf))
    
    @test isnan(xexpy(NaN, -Inf))
    @test isnan(xexpy(NaN, 1))
    @test isnan(xexpy(0, NaN))

    @inferred xexpy(0,-Inf)
    @inferred xexpy(1,1)
    @inferred xexpy(1.,-Inf)

    @test xexpy(2,3) ≈ 2exp(3)

    for x = -10:10, y = -10:10
        @test xexpy(x,y) ≈ x * exp(y)
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
    @test log1mexp(-0.1) ≈ log1p(-exp(-0.1))
    @test log1mexp(-1.0) ≈ log1p(-exp(-1.0))
    @test log1mexp(-10.0) ≈ log1p(-exp(-10.0))
end


@testset "midpoint" begin
    @test iszero(midpoint(-Inf, Inf))
    @test iszero(midpoint(-1., 1.))
    @test iszero(midpoint(nextfloat(-Inf), prevfloat(Inf)))

    @test isnan(midpoint(NaN, 1.))
    @test isnan(midpoint(1., NaN))
    @test isnan(midpoint(NaN, NaN))
    @test isnan(midpoint(1., -1.))
    @test isnan(midpoint(1., -Inf))
    @test isnan(midpoint(Inf, 1.))
    @test isnan(midpoint(Inf, -Inf))

    @test 1. ≤ midpoint(1., 1. + eps()) ≤ 1. + eps()
    @test 1. ≤ midpoint(1., Inf) == prevfloat(Inf) ≤ Inf
    @test -Inf ≤ midpoint(-Inf, 1.) == nextfloat(-Inf) ≤ 1.
end
