using CairoMakie
using Nodiv
using NodivMakie
using Phylo
using Random
using SpatialEcology
using Test

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
