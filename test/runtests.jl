using Iduna
using Documenter: DocMeta, doctest
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
    include("Utils/utils.jl")
    include("IDMapping/id_mapping.jl")
    include("EPLI/EPLI.jl")
    include("ThorAxeMSA/ThorAxeMSA.jl")
    include("MSAExpansion/MSAExpansion.jl")
    include("ResultsValidation/results_validation.jl")
    include("docs_examples.jl")

    # Run doctests
    DocMeta.setdocmeta!(Iduna, :DocTestSetup, :(using Iduna); recursive = true)
    doctest(Iduna)

    # Smoke test for MAPK8 live integration, which is not suitable for GitHub Actions due
    # to potential instability and long runtime.
    if get(ENV, "GITHUB_ACTIONS", "false") == "true"
        @info "Skipping MAPK8 live integration smoke test on GitHub Actions"
    else
        include("Integration/mapk8_smoke.jl")
    end
end

@run_package_tests
