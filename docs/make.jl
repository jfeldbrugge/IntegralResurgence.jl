ENV["JULIA_PROGRESS"] = "false"

using IntegralResurgence
using Documenter

DocMeta.setdocmeta!(IntegralResurgence, :DocTestSetup, :(using IntegralResurgence); recursive=true)

makedocs(;
    modules=[IntegralResurgence],
    authors="Job Feldbrugge <14946916+jfeldbrugge@users.noreply.github.com> and contributors",
    sitename="IntegralResurgence.jl",
    format=Documenter.HTML(;
        canonical="https://jfeldbrugge.github.io/IntegralResurgence.jl",
        edit_link="main",
        assets=String[],
        size_threshold=nothing,
    ),
    pages=[
        "Home" => "index.md",
        "Theory" => "theory.md",
        "Tutorial" => "tutorial.md",
        "Citations" => "citations.md",
        "Reference" => "reference.md",
    ],
)

deploydocs(;
    repo="github.com/jfeldbrugge/IntegralResurgence.jl",
    push_preview = true,
    devbranch="main",
)
