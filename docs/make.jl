using GLV
using Documenter

DocMeta.setdocmeta!(GLV, :DocTestSetup, :(using GLV); recursive=true)

makedocs(
    sitename="GLV.jl",
    modules=[GLV],
    format=Documenter.HTML(),
    pages=["Welcome" => "index.md", "Manual" => "manual.md"],
)

deploydocs(;
    repo="github.com/ismael-lajaaiti/GLV.jl",
)
