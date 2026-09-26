# Plots of Nodiv's `sos_clusters`, the clusters of nodes by the similarity of their SOS
# maps: the clustered similarity heatmap, and the clusters on the tree. Both number and
# colour the clusters of more than one node by the cluster's `labels`, so the two views can
# be read against each other.

"""
    cluster_colors(m)

The colours of `m` clusters, as drawn by [`sos_cluster_heatmap`](@ref) and
[`cluster_tree`](@ref): cycling through the `:tab20` colour map.
"""
function cluster_colors(m)
    c = Makie.to_colormap(:tab20)
    return [c[mod1(i, length(c))] for i in 1:m]
end

"""
    cluster_tree(tree, clusters; title = "", markersize = 12, kwargs...)

The clusters of an `SOSClusters` (from Nodiv's `sos_clusters`) on the tree: a marker at each
node of a cluster of more than one node, one colour per cluster, with a legend. The nodes
on their own are left unmarked, so the co-patterned groups stand out. The cluster numbers
and colours are those outlined on the [`sos_cluster_heatmap`](@ref).

Keyword arguments:
- `treetype = :fan`, `showtips = false`, `strokewidth = 0.5`, `strokecolor = :gray30`,
  and any other keyword arguments: passed to [`treeplot`](@ref)
- `axis = (;)`, `figure = (;)`: attributes for the axis and the figure

Returns a `FigureAxisPlot`.
"""
function cluster_tree(tree, clusters::SOSClusters; title = "", markersize = 12,
                      treetype = :fan, showtips = false, strokewidth = 0.5,
                      strokecolor = :gray30, axis = (;), figure = (;), kwargs...)
    shown = Dict(n => clusters.labels[c] for (n, c) in clusters.groups
                 if haskey(clusters.labels, c))
    fig, ax, p = treeplot(tree; treetype, showtips, nodegroup = shown,
                          groupcolors = cluster_colors(length(clusters.labels)),
                          markersize, strokewidth, strokecolor,
                          axis = merge((; title), axis),
                          figure = merge((; size = (900, 800)), figure), kwargs...)
    isempty(shown) || Legend(fig[1, 2], ax, "cluster")
    return Makie.FigureAxisPlot(fig, ax, p)
end

"""
    sos_cluster_heatmap(clusters; title = "", labelsize = 9, figure = (;)) -> Figure

The clustering of an `SOSClusters` (from Nodiv's `sos_clusters`) as a heatmap of the
similarity `|r|` of the nodes' SOS maps, in the order of the dendrogram, which is drawn on
the left. Genuine groups, if any, show as dense blocks on the diagonal against a uniform
background. The clusters of more than one node are outlined on the diagonal, with their
number and colour from [`cluster_tree`](@ref).

The node names are on both axes; the leaf order, bottom to top, is
`clusters.nodes[clusters.hclust.order]`. For many nodes, make the figure larger or
`labelsize` smaller. `figure` holds attributes for the figure.
"""
function sos_cluster_heatmap(clusters::SOSClusters; title = "", labelsize = 9,
                             figure = (;))
    hc = clusters.hclust
    ord = hc.order
    labs = clusters.nodes[ord]                  # Node names in dendrogram-leaf order
    S = (1 .- clusters.distances)[ord, ord]     # Similarity |r|, reordered to match
    n = length(labs)
    grp = [clusters.groups[node] for node in labs]  # Cluster of each leaf, in leaf order

    fig = Figure(; size = (1300, 1150), figure...)
    Label(fig[0, 1:3], title; fontsize = 18, font = :bold)
    # The dendrogram on the left, its leaves (height 0) against the heatmap; its heights are
    # drawn at negative x, so the ticks show them as positive
    dend = Axis(fig[1, 1]; xlabel = "1 - |r|",
                xtickformat = xs -> string.(round.(abs.(xs); digits = 2)),
                xgridvisible = false, ygridvisible = false, yticksvisible = false,
                yticklabelsvisible = false, leftspinevisible = false,
                topspinevisible = false, rightspinevisible = false)
    hm = Axis(fig[1, 2]; xticks = (1:n, labs), yticks = (1:n, labs),
              xticklabelrotation = pi / 2, xticklabelsize = labelsize,
              yticklabelsize = labelsize)
    # Leaf hc.order[k] at y = k, lined up with heatmap row k; each merge at its height
    dendrogram!(dend, Makie.hcl_nodes(hc; useheight = true); absolute = true,
                rotation = :right, color = :black)
    h = heatmap!(hm, 1:n, 1:n, S; colormap = :viridis, colorrange = (0, 1))
    Colorbar(fig[1, 3], h; label = "|r|")

    colors = cluster_colors(length(clusters.labels))
    for (c, i) in clusters.labels               # Cut clusters are contiguous in `ord`
        rows = findall(==(c), grp)
        lo, hi = extrema(rows)
        box = Rect2d(lo - 0.5, lo - 0.5, hi - lo + 1, hi - lo + 1)
        poly!(hm, box; color = :transparent, strokecolor = :black, strokewidth = 4)
        poly!(hm, box; color = :transparent, strokecolor = colors[i], strokewidth = 2)
        textlabel!(hm, Point2d(hi + 0.5, hi + 0.5); text = string(i), fontsize = 11,
                   background_color = colors[i], strokewidth = 1, padding = 2,
                   cornerradius = 2,
                   # Inside the corner at the top edge, so it is not clipped
                   text_align = hi == n ? (:right, :top) : (:left, :bottom))
    end

    linkyaxes!(dend, hm)
    xlims!(dend, -1.02maximum(hc.heights), 0)
    limits!(hm, 0.5, n + 0.5, 0.5, n + 0.5)
    colsize!(fig.layout, 1, Relative(0.20))
    return fig
end
