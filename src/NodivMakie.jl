module NodivMakie

using Reexport
@reexport using Makie
@reexport using Phylo
@reexport using SpatialEcology
@reexport using Nodiv
using Makie: automatic, Automatic, Colorant, RGBAf, Point2d, Vec2f, to_color
using ColorTypes: alpha
import EcoBase
using SparseArrays: SparseMatrixCSC, sparse, rowvals, nonzeros, nzrange
using EcoBase: richness, places, getcoords, coordinates, xrange, yrange, occupancy
import FileIO, ImageIO
using Makie: ncols, insertcols!
import MultivariateStats

include("layout.jl")
include("treeplot.jl")
include("sitemap.jl")
include("nodepanel.jl")
include("images.jl")
include("imageclades.jl")
include("treeimages.jl")
include("ordination.jl")
include("explorer.jl")

export TreeLayout, treelayout
export treeplot, treeplot!, TreePlot
export sitemap, sitemap!, SiteMap
export nodepanel, nodepanel!, NodePanel, hassos, cladecolors
export nodeexplorer, NodeExplorer, nodeat, onnodeclick, focuscolors
export SpeciesImages, speciesname, treeimages!, TreeImages, missingimages, cladeimages!
export imagegeometry, imageclades, selectclades, CladeImage
export sosordination, SOSOrdination, ordinationplot, ordinationplot!, OrdinationPlot

end
