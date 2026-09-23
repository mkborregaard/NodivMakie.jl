module NodivMakie

using Reexport
@reexport using Makie
@reexport using Phylo
@reexport using SpatialEcology
@reexport using Nodiv
using Makie: automatic, Automatic, Colorant, RGBAf, Point2d, Vec2f, to_color

include("layout.jl")
include("treeplot.jl")

export TreeLayout, treelayout
export treeplot, treeplot!, TreePlot

end
