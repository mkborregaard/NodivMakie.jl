module NodivMakie

using ColorTypes: alpha
using EcoBase: EcoBase
using EcoBase: coordinates, getcoords, occupancy, places, richness, xrange, yrange
using FileIO: FileIO
using GridLayoutBase: GridLayoutBase
using GridLayoutBase: insertcols!
using ImageIO: ImageIO      # loaded here so that FileIO can read PNG and JPEG images
using Makie: Automatic, Colorant, Point2d, RGBAf, Vec2f
using Makie: automatic, ncols, to_color
using Reexport: @reexport

@reexport using Makie
@reexport using Nodiv
@reexport using Phylo
@reexport using SpatialEcology

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
include("compat.jl")

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
