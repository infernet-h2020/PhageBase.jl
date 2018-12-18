module PhageBase

using PhageSeq

export Sequence, SequencesIterator, hamming, subseq, seqinsert # from PhageSeq
export energy

include("util.jl")
include("math.jl")

include("fermi_dirac.jl")

include("FastSeq.jl")

include("Fields.jl")
include("FieldsChem.jl")

include("prior.jl")

include("Dataset.jl")

include("gauges.jl")

include("clean.jl")

end # module
