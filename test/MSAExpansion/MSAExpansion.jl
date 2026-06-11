@testset "MSAExpansion" begin
    mktempdir() do _tmp
        global tmp = _tmp
        include("stockholm_and_alignment.jl")
        include("cache_and_state.jl")
        include("workflow.jl")
    end
end
