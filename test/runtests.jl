using NodivMakie
using CairoMakie
using Test
using Random

include("testdata.jl")

@testset "NodivMakie.jl" begin
    include("treeplot.jl")
    include("sitemap.jl")
    include("nodepanel.jl")
    include("ordination.jl")
    include("explorer.jl")
    include("figures.jl")
    include("clusters.jl")
    include("images.jl")
end
