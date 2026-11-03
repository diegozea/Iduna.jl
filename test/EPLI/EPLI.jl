@testset "EPLI" begin
    include("fixtures.jl")
    include("score_types.jl")
    include("core_scoring.jl")
    include("helpers_and_aligners.jl")
    include("normalization_arithmetic.jl")
    include("jll_extensions.jl")
end
