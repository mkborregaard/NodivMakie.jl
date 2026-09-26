module NodivMakie

using ColorTypes: Colorant, alpha
using EcoBase: EcoBase
using EcoBase: coordinates, getcoords, occupancy, places, richness, xrange, yrange
using FileIO: FileIO, save
using GridLayoutBase: GridLayoutBase, insertcols!
using ImageIO: ImageIO      # loaded here so that FileIO can read PNG and JPEG images
using Makie: Makie
using Makie: Axis, Colorbar, Consume, DataInspector, Figure, GridLayout, Label, Legend
using Makie: Lines, Mouse, Observable, Point2d, Point2f, RGBAf, Rect2d, Relative, Reverse
using Makie: Scatter, Vec2f
using Makie: @recipe
using Makie: (..), automatic, barplot!, colsize!, dendrogram!, events, heatmap!
using Makie: hidedecorations!, hidespines!, image!, lift, limits!, lines!, linkaxes!
using Makie: linkyaxes!, mouseposition, mouseposition_px, ncols, on, pick, poly!
using Makie: scatter!, text!, textlabel!, to_color, to_colormap, to_value, translate!
using Makie: xlims!
using Nodiv: Nodiv
using Nodiv: SOSClusters, SOSOrdination
using Nodiv: clade_richness, default_score, divergent_nodes, most_divergent, node_scores
using Nodiv: nodespecies, sos_ordination
using Phylo: Phylo
using Phylo: getchildren, getdescendants, getnodedata, getnodename, getparent, hasinbound
using Phylo: hasnode, isleaf
using SpatialEcology: sitestats

export TreeLayout, tree_layout
export treeplot, treeplot!, TreePlot
export sitemap, sitemap!, SiteMap
export node_panel, node_panel!, NodePanel, has_sos, clade_colors, sos_map!
export node_explorer, NodeExplorer, node_at, on_node_click, focus_colors, link_explorers!
export explorer_tree!, ExplorerTree
export SpeciesImages, species_name, tree_images!, TreeImages, missing_images, clade_images!
export image_geometry, image_clades, select_clades, CladeImage
export ordinationplot, ordinationplot!, OrdinationPlot
export eigenvalueplot, eigenvalueplot!, EigenvaluePlot
export map_figure, metric_tree, node_panel_pdf
export cluster_colors, cluster_tree, sos_cluster_heatmap

# Deprecated
export cladecolors, cladeimages!, explorertree!, focuscolors, hassos, imageclades
export imagegeometry, missingimages, nodeat, nodeexplorer, nodepanel, nodepanel!
export onnodeclick, selectclades, sosmap!, speciesname, treeimages!, treelayout

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
include("deprecated.jl")

end
