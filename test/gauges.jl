Random.seed!(1)
A = 4; L = 3;
flen = fieldslen(A,L)
fieldschem = FieldsChem{A,L}(randn(flen + 1))
f1 = zero_sum_gauge(fieldschem)

@test length(fieldschem) == length(f1) == flen + 1
@test is_zero_sum_gauge(f1; atol=1e-10)

for s in SequencesIterator{A,L}()
    @test bindφ(fieldschem, s) ≈ bindφ(f1, s)
end


# a longer sequence

A = 21; L = 35;
flen = fieldslen(A,L)
fieldschem = FieldsChem{A,L}(randn(flen + 1))
f1 = zero_sum_gauge(fieldschem)

@test length(fieldschem) == length(f1) == flen + 1
@test is_zero_sum_gauge(f1; atol=1e-10)

for _ = 1:20
    s = rand(Sequence{A,L})
    @test bindφ(fieldschem, s) ≈ bindφ(f1, s)
end
