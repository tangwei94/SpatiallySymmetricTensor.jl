using IPEPSC6v
using Documenter

DocMeta.setdocmeta!(IPEPSC6v, :DocTestSetup, :(using IPEPSC6v); recursive=true)

makedocs(;
    modules=[IPEPSC6v],
    authors="Wei Tang <tangwei@smail.nju.edu.cn> and contributors",
    repo="https://github.com/tangwei94/IPEPSC6v.jl/blob/{commit}{path}#{line}",
    sitename="IPEPSC6v.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://tangwei94.github.io/IPEPSC6v.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/tangwei94/IPEPSC6v.jl",
    devbranch="main",
)
