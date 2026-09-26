module NodivMakie

using Reexport
@reexport using Makie
@reexport using Phylo
@reexport using SpatialEcology
@reexport using Nodiv
using Makie: automatic, Automatic, Colorant, RGBAf, Point2d, Vec2f, to_color
using ColorTypes: alpha
import EcoBase
using EcoBase: richness, places, getcoords, coordinates, xrange, yrange, occupancy
import FileIO, ImageIO
using Makie: ncols, insertcols!

include("layout.jl")
include("treeplot.jl")
include("sitemap.jl")
include("nodepanel.jl")
include("images.jl")
include("imageclades.jl")
include("treeimages.jl")
include("ordination.jl")
include("explorer.jl")
include("figures.jl")
include("clusters.jl")

export TreeLayout, treelayout
export treeplot, treeplot!, TreePlot
export sitemap, sitemap!, SiteMap
export nodepanel, nodepanel!, NodePanel, hassos, cladecolors, sosmap!
export nodeexplorer, NodeExplorer, nodeat, onnodeclick, focuscolors, link_explorers!
export explorertree!, ExplorerTree
export SpeciesImages, speciesname, treeimages!, TreeImages, missingimages, cladeimages!
export imagegeometry, imageclades, selectclades, CladeImage
export ordinationplot, ordinationplot!, OrdinationPlot
export eigenvalueplot, eigenvalueplot!, EigenvaluePlot
export map_figure, metric_tree, node_panel_pdf
export cluster_colors, cluster_tree, sos_cluster_heatmap

end
