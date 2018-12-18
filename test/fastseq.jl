@testset "FastSeq" begin
    for testrep=1:3
        Random.seed!(476272059+testrep)

        A = rand(2:6); L = rand(2:6)
        s = rand(Sequence{A,L})

        @test FastSeq{3}(1,2) == FastSeq{4}(1,2)
        @test FastSeq{3}(1,2) != FastSeq{4}(1,3)

        @test FastSeq{3}(1,2) == Sequence{4}(1,2)
        @test FastSeq{3}(1,2) != Sequence{4}(1,3)

        fs = FastSeq(s)
        @test fs.sequence == s
        @test length(fs.fieldidx) == binomial(L+1,2)
        @test issorted(fs.fieldidx)

        for i=1:L
            @test (i-1)A+s[i] ∈ fs.fieldidx
        end

        for i=1:L, j=i+1:L
            @test A*L + binomial(j-1,2)*A^2 + (i-1)A^2 + (s[j]-1)A + s[i] ∈ fs.fieldidx
        end
    end
end
