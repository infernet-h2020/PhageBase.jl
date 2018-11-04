@testset "train split" begin
    A = 5; L = 5; V = 2; T = 3; N0 = 100
    f = FieldsChem{A,L}()
    randn!(f.x)
    seqs = collect(SequencesIterator{A,L}())
    S = length(seqs)
    N = rand(1:10^6, S, V, T)
    data = Dataset(seqs, N)

    traindata, testsdata = randtrain(data, 1000)

    @test sort(union(traindata.sequences, testsdata.sequences)) == sort(data.sequences)
    @test isempty(intersect(traindata.sequences, testsdata.sequences))
end