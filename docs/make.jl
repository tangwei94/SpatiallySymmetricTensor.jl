using SpatiallySymmetricTensors
using Documenter

DocMeta.setdocmeta!(SpatiallySymmetricTensors, :DocTestSetup, :(using SpatiallySymmetricTensors); recursive=true)

makedocs(;
    modules=[SpatiallySymmetricTensors],
    authors="Wei Tang <tangwei@smail.nju.edu.cn> and contributors",
    repo="https://github.com/tangwei94/SpatiallySymmetricTensors.jl/blob/{commit}{path}#{line}",
    sitename="SpatiallySymmetricTensors.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://tangwei94.github.io/SpatiallySymmetricTensors.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/tangwei94/SpatiallySymmetricTensors.jl",
    devbranch="main",
)
