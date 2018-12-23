module PhageBase

using Random, ArgCheck

export Sequence, SequencesIterator, hamming, subseq, seqinsert
export energy

include("util.jl")
include("math.jl")

include("Sequence.jl")
include("SequencesIterator.jl")

include("fermi_dirac.jl")

include("FastSeq.jl")

include("Fields.jl")
include("FieldsChem.jl")

include("prior.jl")

include("Dataset.jl")

include("gauges.jl")

include("clean.jl")

end # module
