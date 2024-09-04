using GLV
using Documenter
using Literate

# Literate.markdown(
#     joinpath(@__DIR__, "src", "cavity.jl"), joinpath(@__DIR__, "src");
#     credit=false
# )

DocMeta.setdocmeta!(GLV, :DocTestSetup, :(using GLV); recursive=true)

makedocs(
    sitename="GLV.jl",
    modules=[GLV],
    format=Documenter.HTML(
        prettyurls=get(ENV, "CI", nothing) == "true"
    ),
    pages=[
        "Home" => "index.md",
        "The cavity method" => "cavity.md",
    ],
)

deploydocs(;
    repo="github.com/ismael-lajaaiti/GLV.jl",
)
