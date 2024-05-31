using SpatiallySymmetricTensor
using Documenter

DocMeta.setdocmeta!(SpatiallySymmetricTensor, :DocTestSetup, :(using SpatiallySymmetricTensor); recursive=true)

makedocs(;
    modules=[SpatiallySymmetricTensor],
    authors="Wei Tang <tangwei@smail.nju.edu.cn> and contributors",
    repo="https://github.com/tangwei94/SpatiallySymmetricTensor.jl/blob/{commit}{path}#{line}",
    sitename="SpatiallySymmetricTensor.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://tangwei94.github.io/SpatiallySymmetricTensor.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/tangwei94/SpatiallySymmetricTensor.jl",
    devbranch="main",
)
