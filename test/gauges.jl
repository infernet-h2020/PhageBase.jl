@testset "gauges" begin
    Random.seed!(1)
    A = 4; L = 3;
    flen = fieldslen(A,L)
    f = FieldsChem{A,L}(randn(flen + 1))
    f1 = zero_sum_gauge(f)

    @test length(f) == length(f1) == flen + 1
    @test is_zero_sum_gauge(f1; atol=1e-10)
    for s in SequencesIterator{A,L}()
        @test bindφ(f, s) ≈ bindφ(f1, s)
    end

    
    A = 21; L = 35;
    flen = fieldslen(A,L)
    f = FieldsChem{A,L}(randn(flen + 1))
    f1 = zero_sum_gauge(f)

    @test length(f) == length(f1) == flen + 1
    @test is_zero_sum_gauge(f1; atol=1e-10)

    for _ = 1:20
        s = rand(Sequence{A,L})
        @test bindφ(f, s) ≈ bindφ(f1, s)
    end
end