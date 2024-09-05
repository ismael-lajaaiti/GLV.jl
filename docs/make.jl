using GLV
using Documenter
using DocumenterCitations
using Literate

# Generate markdown files from julia scripts.
Literate.markdown(
    joinpath(@__DIR__, "src", "cavity.jl"), joinpath(@__DIR__, "src");
    credit=false
)

DocMeta.setdocmeta!(GLV, :DocTestSetup, :(using GLV); recursive=true)

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style=:authoryear
)

makedocs(
    sitename="GLV.jl",
    modules=[GLV],
    format=Documenter.HTML(
        prettyurls=get(ENV, "CI", nothing) == "true"
    ),
    pages=[
        "Home" => "index.md",
        "The cavity method" => "cavity.md",
        "References" => "references.md",
    ],
    plugins=[bib],
)

deploydocs(;
    repo="github.com/ismael-lajaaiti/GLV.jl",
)
