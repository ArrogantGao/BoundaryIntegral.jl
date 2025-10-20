using BoundaryIntegral
using Documenter

DocMeta.setdocmeta!(BoundaryIntegral, :DocTestSetup, :(using BoundaryIntegral); recursive=true)

makedocs(;
    modules=[BoundaryIntegral],
    authors="Xuanzhao Gao <xgao@flatironinstitute.org> and contributors",
    sitename="BoundaryIntegral.jl",
    format=Documenter.HTML(;
        canonical="https://Xuanzhao Gao.github.io/BoundaryIntegral.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Xuanzhao Gao/BoundaryIntegral.jl",
    devbranch="main",
)
