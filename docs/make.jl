using GLV
using Documenter
using DocumenterCitations
using Literate

# Generate markdown files from julia scripts.
function generate_md(filename)
    Literate.markdown(
        joinpath(@__DIR__, "src", filename),
        joinpath(@__DIR__, "src");
        credit = false,
    )
end
generate_md("cavity.jl")
generate_md("functional-extinctions.jl")
generate_md("response-to-pulse.jl")

DocMeta.setdocmeta!(GLV, :DocTestSetup, :(using GLV); recursive = true)

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style = :authoryear)

makedocs(;
    sitename = "GLV.jl",
    modules = [GLV],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = [
        "Home" => "index.md",
        "The cavity method" => "cavity.md",
        "Interactions structure species responses" => "response-to-pulse.md",
        "Functional extinctions" => "functional-extinctions.md",
        "References" => "references.md",
    ],
    plugins = [bib],
)

deploydocs(; repo = "github.com/ismael-lajaaiti/GLV.jl")
