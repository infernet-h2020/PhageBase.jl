using Random

@testset "fields, energy" begin
    for testrep=1:5
        Random.seed!(686371329+testrep)

        A = rand(2:6); L = rand(2:6); 
        V = rand(2:6); T = rand(2:6);

        h = Dict((a,i) => randn() for a = 1:A for i = 1:L)
        J = Dict((a,b,i,j) => randn() for a = 1:A for b = 1:A for i = 1:L for j = i+1:L)
        
        fieldsh = [h[(a,i)] for i = 1:L for a = 1:A]
        fieldsJ = [J[(a,b,i,j)] for j = 2:L for i = 1:j-1 for b = 1:A for a = 1:A]
        
        @test PhageBase.fieldslen(A,L) == length(fieldsh) + length(fieldsJ)

        fields = Fields{A,L}([fieldsh; fieldsJ])
        
        for i = 1:L, a = 1:A
            @test fields[a,i] == h[(a,i)]
        end
        for i = 1:L, j = i+1 : L, a = 1:A, b = 1:A
            @test fields[a,b,i,j] == J[(a,b,i,j)]
        end

        s = rand(Sequence{A,L})
        Hh = -sum(h[(s[i], i)] for i = 1:L)
        HJ = -sum(J[(s[i], s[j], i, j)] for i = 1:L for j = i+1:L)

        @test energy(fields, s) ≈ Hh + HJ
        @test energy(fields, FastSeq(s)) ≈ Hh + HJ


        @test (fields[2] = 1.) == 1.
        @test fields[2] == 1.

        @test (fields[2,1] = 2.) == 2.
        @test fields[2,1] == 2.

        @test (fields[2,1,1,2] = 3.) == 3.
        @test fields[2,1,1,2] == 3.
    end
end


@testset "FieldsChem" begin
    Random.seed!(42572)
    A = rand(2:6); L = rand(2:6);
    flen = fieldslen(A,L)
    fields = FieldsChem{A,L}(randn(flen + 1))
    @test fields.μ == fields[end] == fields.x[end]
    @test length(fields) == flen + 1

    oldμ = fields.μ
    fields.μ += 1
    @test fields.μ == oldμ + 1

    f0 = Fields(fields)
    @test f0.x === fields.x

    for _ = 1:5
        s = rand(Sequence{A,L})
        @test energy(fields, s) == energy(f0, s)
        @test bindφ(fields, s) == fields.μ - energy(fields, s)
        @test fermi_dirac_prob(fields, s) == fermi_dirac_prob(fields.μ, energy(fields, s))
    end
end