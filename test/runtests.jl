using Iduna
using Test
using TestItemRunner

@testset "Iduna.jl" begin
    include("aqua.jl")
    include("complexity.jl")
    include("namespace.jl")
    include("App/app.jl")
    include("API/API.jl")
    include("example_fixtures.jl")
    include("readme_examples.jl")
    include("docs_examples.jl")
    include("Utils/utils.jl")
    include("IDMapping/id_mapping.jl")
    include("EPLI/EPLI.jl")
    include("ThorAxeMSA/ThorAxeMSA.jl")
    include("MSAExpansion/MSAExpansion.jl")
    include("ResultsValidation/results_validation.jl")
    if get(ENV, "GITHUB_ACTIONS", "false") == "true"
        @info "Skipping MAPK8 live integration smoke test on GitHub Actions"
    else
        include("Integration/mapk8_smoke.jl")
    end
end

@run_package_tests
