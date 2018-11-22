
Random.seed!(42572)

A = rand(2:6); L = rand(2:6);
flen = fieldslen(A,L)
fields = FieldsChem{A,L}(randn(flen + 1))

@test getμ(fields) == fields.x[end]
@test length(fields) == flen + 1

setμ!(fields, 1e5); @test getμ(fields) == 1e5

s = rand(Sequence{A,L})
@test energy(fields, s) == energy(Fields(fields), s)
@test bindφ(fields, s) == getμ(fields) - energy(fields, s)
@test fermi_dirac_prob(fields, s) ==
    fermi_dirac_prob(getμ(fields), energy(fields, s)) ==
    fermi_dirac_prob(bindφ(fields, s))
