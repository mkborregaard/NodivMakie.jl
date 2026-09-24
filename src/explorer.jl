# Linking a tree to other figure elements: picking a node on a `treeplot`, and the
# node explorer (tree on the left, node panel on the right, click to select).

"""
    nodeat(treeplot, plot, index)

The node name under a pick result `(plot, index)` of `pick(scene, ...)`, or `nothing`.
A node marker gives its node, and a branch gives the node it leads to.
"""
function nodeat(tp::TreePlot, plt, idx)
    (plt === nothing || idx === nothing || idx < 1) && return nothing
    plt.parent === tp || return nothing
    names = tp.tree_layout[].names
    if plt isa Lines
        idx <= length(tp.branch_owner[]) || return nothing
        return names[tp.branch_owner[][idx]]
    elseif plt isa Scatter
        pts = plt[1][]
        idx <= length(pts) || return nothing
        # markers (plain or grouped) sit on node positions; the padding points do not
        i = findfirst(==(pts[idx]), tp.node_points[])
        return i === nothing ? nothing : names[i]
    end
    return nothing
end

"""
    onnodeclick(f, ax, treeplot; range = 10)

Call `f(nodename)` when a node marker or branch of `treeplot` in `ax` is left-clicked.
`range` is the pick radius in pixels. Needs an interactive backend (GLMakie, WGLMakie).
"""
onnodeclick(f, ax::Axis, tp::TreePlot; range = 10) = _onnodeclick(f, ax, tp, pick; range)

# `pickfn` is the backend's `pick`; the tests swap in a stand-in, as CairoMakie cannot pick
function _onnodeclick(f, ax::Axis, tp::TreePlot, pickfn; range = 10)
    scene = ax.scene
    return on(events(scene).mousebutton; priority = 2) do ev
        (ev.button == Mouse.left && ev.action == Mouse.press) || return Consume(false)
        Makie.is_mouseinside(scene) || return Consume(false)
        n = nodeat(tp, pickfn(scene, events(scene).mouseposition[], range)...)
        n === nothing && return Consume(false)
        f(n)
        return Consume(true)
    end
end

"""
    focuscolors(tree, layout, node, sos_colormap, contextcolor; inset = 0.15)

One branch colour per node of `layout` (a [`TreeLayout`](@ref)), for `treeplot`'s
`branchcolor`. The clade of `node`'s first child is drawn in the high end of
`sos_colormap`, the clade of its second child in the low end, and every other branch
in `contextcolor`. The two colours are taken `inset` in from the ends of the colour map
(0 = the end colours), which for a diverging map like RdYlBu gives lighter colours.

This follows Nodiv's SOS: a node's SOS is the standardised richness of its FIRST child
(`getchildren(tree, node)[1]`) against the null, which keeps the parent clade's richness
per cell fixed. Positive SOS (the high end) marks cells where the first child is
over-represented; negative SOS (the low end) marks cells where the second child is.
"""
function focuscolors(tree, layout::TreeLayout, node, sos_colormap, contextcolor;
                     inset = 0.15)
    cmap = Makie.to_colormap(sos_colormap)
    high = Makie.interpolated_getindex(cmap, 1.0 - inset)
    low = Makie.interpolated_getindex(cmap, Float64(inset))
    colors = fill(to_color(contextcolor), length(layout))
    for (child, color) in zip(getchildren(tree, node)[1:2], (high, low))
        colors[layout.index[getnodename(tree, child)]] = color
        for d in getdescendants(tree, child)
            colors[layout.index[getnodename(tree, d)]] = color
        end
    end
    return colors
end

"""
    NodeExplorer

The parts of a [`nodeexplorer`](@ref): `figure`, the tree `axis` and `treeplot`, the
[`NodePanel`](@ref) `panel`, `status`, the text of the label above the tree, and
`images`, the species images around the tree ([`TreeImages`](@ref)) or `nothing`.
"""
struct NodeExplorer
    figure::Figure
    axis::Axis
    treeplot::TreePlot
    panel::NodePanel
    status::Observable{String}
    images::Union{TreeImages, Nothing}
end

metricvalues(res, metric::Symbol) = getfield(res, metric)
metricvalues(res, metric::AbstractDict) = metric

# The metric shown by default: Nodiv's default divergence score for the result type
defaultmetric(::Nodiv.NodeMetrics) = :rms
defaultmetric(res) = :gnd

# The nodes marked by default: Nodiv's `divergent_nodes`, with its default threshold,
# ranked by the metric shown where `divergent_nodes` supports it
function defaultnodes(res::Nodiv.NodeMetrics, metric)
    by = metric isa Symbol && metric in (:rms, :pval, :gnd) ? metric : :rms
    return divergent_nodes(res; by)
end
defaultnodes(res::Nodiv.NodeAnalysis, metric) = divergent_nodes(res)
defaultnodes(res, metric) = collect(keys(sosvalues(res)))   # a plain SOS Dict

"""
    nodeexplorer(assemblage, tree, res; kwargs...)

The tree and a [`nodepanel`](@ref) side by side. The tree marks `nodes` coloured by
`metric`; clicking a marker or a branch shows that node in the panel. The two clades
below the node shown are drawn in colours from the two ends of the SOS colour map (see
[`focuscolors`](@ref)), and the rest of the tree in `contextcolor`. `res` is the cached
`NodeMetrics`/`NodeAnalysis`.

Keyword arguments:
- `metric`: a field of `res` (`:rms`, `:gnd`, `:pval`, `:ses`, ...) or a Dict of node =>
  value. Default `:rms` for a `NodeMetrics`, `:gnd` for a `NodeAnalysis` (Nodiv's defaults).
- `nodes`: the nodes to mark. Default: the divergent nodes, `divergent_nodes(res)` with
  Nodiv's default threshold, ranked by `metric` when that is `:rms`, `:pval` or `:gnd`.
  Pass `divergent_nodes(res; threshold = ...)` for another threshold, or `:all` for every
  node with an SOS.
- `node`: the node shown first; default the most divergent marked node (highest metric
  value, lowest for `:pval`)
- `treetype = :fan`, `showtips = false`, `colormap = :YlOrRd`, `colorrange = automatic`,
  `markersize = 14`, `strokewidth = 1`, `strokecolor = :gray20` (the marker outline):
  passed to `treeplot`, as is everything in `treekw = (;)`. The branch colours are set by
  the explorer, so `branchcolor` is overridden.
- `contextcolor = :gray75`: the branches outside the two clades of the node shown
- `focusinset = 0.15`: how far in from the ends of the SOS colour map the two clade
  colours are taken; 0 gives the end colours, larger is lighter
- `images = nothing`: species images to draw around the tree, as a
  [`SpeciesImages`](@ref) or a directory path; see [`treeimages!`](@ref). Range sizes
  are taken from `assemblage` unless `imageoptions` has a `rangesize`.
- `imageoptions = (;)`: keyword arguments for `treeimages!`
- `panel = (;)`: keyword arguments for `nodepanel!`
- `figure = (;)`: attributes for the `Figure`
- `pickfn = pick`: the picking function (replaced in the tests, where CairoMakie cannot pick)

Returns `(figure, explorer)`; see [`NodeExplorer`](@ref). Clicking needs an interactive
backend (`using GLMakie`); with CairoMakie you get the figure for the initial node.
"""
function nodeexplorer(assemblage, tree, res; metric = defaultmetric(res),
                      nodes = automatic, node = automatic, treetype = :fan,
                      showtips = false, colormap = :YlOrRd, colorrange = automatic,
                      markersize = 14, strokewidth = 1, strokecolor = :gray20,
                      contextcolor = :gray75, focusinset = 0.15, images = nothing,
                      imageoptions = (;), treekw = (;), panel = (;), figure = (;),
                      pickfn = pick)
    sos = sosvalues(res)
    vals = metricvalues(res, metric)
    nodes = nodes === automatic ? defaultnodes(res, metric) :
            nodes === :all ? collect(keys(sos)) : collect(nodes)
    marked = Dict(n => Float64(vals[n]) for n in nodes
                  if haskey(vals, n) && isfinite(vals[n]) && hassos(tree, sos, n))
    isempty(marked) && throw(ArgumentError("None of `nodes` has both an SOS and a finite " *
                                           "metric value; pass `nodes = :all` to mark every node with an SOS"))
    mostdivergent = metric === :pval ? argmin : argmax
    node = node === automatic ? mostdivergent(n -> marked[n], keys(marked)) : String(node)
    label = metric isa Symbol ? string(metric) : "value"

    fig = Figure(; size = (1600, 850), figure...)
    treegrid = fig[1, 1] = GridLayout()
    status = Observable("")
    Label(treegrid[1, 1], status; tellwidth = false, halign = :left)
    ax = Axis(treegrid[2, 1]; autolimitaspect = treetype === :fan ? 1 : nothing)
    hidedecorations!(ax)
    hidespines!(ax)
    np = nodepanel!(fig[1, 2], assemblage, tree, node, res; panel...)
    colsize!(fig.layout, 1, Relative(0.45))

    # the node shown: its two clades in the SOS colours, the rest of the tree greyed
    layout = treelayout(tree)
    sos_colormap = np.maps[2].colormap[]
    focus(n) = focuscolors(tree, layout, n, sos_colormap, contextcolor; inset = focusinset)
    # markers only on the marked nodes: a NaN colour hides a marker's fill, not its outline
    tp = treeplot!(ax, tree; treetype, showtips, nodecolor = marked,
                   shownodes = collect(keys(marked)), markersize,
                   strokewidth, strokecolor, colormap, colorrange,
                   merge(treekw, (; branchcolor = focus(node)))...)
    on(n -> (tp.branchcolor = focus(n)), np.node)
    ti = nothing
    if images !== nothing
        rangesize = get(imageoptions, :rangesize, assemblage)
        opts = Base.structdiff(imageoptions, NamedTuple{(:rangesize,)})
        ti = treeimages!(ax, tp, images, rangesize; opts...)
    end
    Colorbar(treegrid[3, 1], tp; vertical = false, flipaxis = false, label,
             tellheight = true, width = Relative(0.6))

    describe(n) = haskey(vals, n) && vals[n] isa Real ? "$n   $label = $(round(vals[n]; sigdigits = 3))" : n
    status[] = describe(node)
    on(n -> (status[] = describe(n)), np.node)
    _onnodeclick(ax, tp, pickfn) do n
        if hassos(tree, sos, n)
            np.node[] = n
        else
            status[] = "$n: no SOS (a tip, or not analysed)"
        end
    end
    return fig, NodeExplorer(fig, ax, tp, np, status, ti)
end
