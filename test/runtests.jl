using PhageBase
using Test, Random
using ForwardDiff

@testset "Fields" begin include("Fields.jl") end
@testset "FieldsChem" begin include("FieldsChem.jl") end
@testset "math" begin include("math.jl") end
@testset "fermi_dirac" begin include("fermi_dirac.jl") end
@testset "sequences" begin include("sequence.jl") end
@testset "prior" begin include("prior.jl") end
@testset "gauges" begin include("gauges.jl") end
@testset "data" begin include("data.jl") end
