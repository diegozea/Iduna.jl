import JSON
import Logging

@testset "API" begin
    include("core_and_fixtures.jl")
    include("load_result.jl")
    include("stage_helpers_and_logging.jl")
    include("pretty_printing.jl")
    include("paths_and_options.jl")
    include("expansion_modes.jl")
    include("failure_artifacts.jl")
end
