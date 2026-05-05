using Iduna
using Test
using TestItemRunner

@testset "Iduna.jl" begin
    include("aqua.jl")
    include("App/app.jl")
    include("API/api.jl")
    include("Utils/utils.jl")
    include("IDMapping/id_mapping.jl")
    include("ThorAxeMSA/thoraxe_msa.jl")
    include("MSAExpansion/msa_expansion.jl")
    include("ResultsValidation/results_validation.jl")
end

@run_package_tests
