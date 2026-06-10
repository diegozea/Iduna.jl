using Pkg

# The docs manifest is ignored, so make the documented package available when
# the docs command is run from a fresh checkout.
cd(@__DIR__) do
    Pkg.develop(PackageSpec(path = ".."))
    Pkg.instantiate()
end

using Iduna
using Documenter

DocMeta.setdocmeta!(Iduna, :DocTestSetup, :(using Iduna); recursive = true)

makedocs(;
    modules = [Iduna, Iduna.IDMapping, Iduna.EPLI, Iduna.ThorAxeMSA,
        Iduna.MSAExpansion, Iduna.ResultsValidation, Iduna.Utils],
    authors = "Diego Javier Zea <diegozea@gmail.com> and contributors",
    sitename = "Iduna.jl",
    format = Documenter.HTML(;
        canonical = "https://diegozea.github.io/Iduna.jl",
        edit_link = "main",
        assets = String[]
    ),
    pages = [
        "Home" => "index.md",
        "EPLI Score" => "epli.md",
        "API" => "api.md",
        "Output Layout" => "output.md"
    ]
)

deploydocs(;
    repo = "github.com/diegozea/Iduna.jl.git",
    devbranch = "main"
)
