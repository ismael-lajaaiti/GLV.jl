using GLV
using Documenter

makedocs(
    sitename="GLV.jl",
    modules=[GLV],
    format=Documenter.HTML()
)

deploydocs(;
    repo="github.com/ismael-lajaaiti/GLV.jl",
)
