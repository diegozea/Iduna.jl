using Iduna
using Documenter

DocMeta.setdocmeta!(Iduna, :DocTestSetup, :(using Iduna); recursive=true)

makedocs(;
    modules=[Iduna],
    authors="Diego Javier Zea <diegozea@gmail.com> and contributors",
    sitename="Iduna.jl",
    format=Documenter.HTML(;
        canonical="https://diegozea.github.io/Iduna.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/diegozea/Iduna.jl",
    devbranch="main",
)
