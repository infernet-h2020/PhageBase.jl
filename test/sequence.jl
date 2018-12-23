@test eltype(Sequence{23,4}) == Int

s = Sequence{10,10}(collect(1:10))
for (i,a) in enumerate(s)
    @test a == i
end

@test Sequence{3}(1,2) == Sequence{4}(1,2)
@test Sequence{3}(1,2) != Sequence{4}(1,3)

@test hash(Sequence{3}(1,2)) == hash(Sequence{4}(1,2))

@test collect(SequencesIterator{2,3}()) == [Sequence{2}(1,1,1), Sequence{2}(1,1,2),
                                            Sequence{2}(1,2,1), Sequence{2}(1,2,2),
                                            Sequence{2}(2,1,1), Sequence{2}(2,1,2),
                                            Sequence{2}(2,2,1), Sequence{2}(2,2,2)]

@test length(collect(SequencesIterator{4,5}())) == 4^5
@test length(SequencesIterator{4,5}()) == 4^5

@test hamming(Sequence{4,2}(1,2), Sequence{4,2}(1,2)) == 0
@test hamming(Sequence{4,2}(1,2), Sequence{5,2}(2,3)) == 2

@test collect(Sequence{3,3}(1,2,3)) == [1,2,3]

@test Sequence{3}(1,2,3) ≤ Sequence{3}(2,2,3)
@test Sequence{3}(1,2,2) < Sequence{3}(3,3,3)

@test subseq(Sequence{3}(1,2,3), 1, 2) == Sequence{3}(1,2)
@test seqinsert(Sequence{3}(1,2,3), 1, Sequence{3}(3,3)) == Sequence{3}(3,3,3)
