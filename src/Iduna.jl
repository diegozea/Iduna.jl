module Iduna

export IDMapping,
       MSAExpansion,
       ResultsValidation,
       ThorAxeMSA,
       Utils,
       ExpansionResult,
       IdunaResult,
       ResolvedTarget,
       SeedSelection,
       ThorAxeMSAResult,
       ValidationResult,
       iduna,
       load_expanded_msa,
       load_seed_msa

include("Utils/Utils.jl")
include("IDMapping/IDMapping.jl")
include("ThorAxeMSA/ThorAxeMSA.jl")
include("MSAExpansion/MSAExpansion.jl")
include("ResultsValidation/ResultsValidation.jl")

using .Utils: ExpansionResult, IdunaResult, ResolvedTarget, SeedSelection, ThorAxeMSAResult,
              ValidationResult

include("API.jl")
include("App.jl")

end
