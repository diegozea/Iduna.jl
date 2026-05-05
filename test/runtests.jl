using Iduna
using Test
using TestItemRunner

@testset "Iduna.jl" begin
    include("aqua.jl")
    include("App/runtests.jl")
    include("API/runtests.jl")
    include("Utils/runtests.jl")
    include("IDMapping/runtests.jl")
    include("ThorAxeMSA/runtests.jl")
    include("MSAExpansion/runtests.jl")
    include("ResultsValidation/runtests.jl")
end

@run_package_tests
