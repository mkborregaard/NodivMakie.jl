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
    focuscolors(tree, layout, node, sos_colormap, contextcolor)

One branch colour per node of `layout` (a [`TreeLayout`](@ref)), for `treeplot`'s
`branchcolor`. The clade of `node`'s first child is drawn in the high end of
`sos_colormap`, the clade of its second child in the low end, and every other branch
in `contextcolor`.

This follows Nodiv's SOS: a node's SOS is the standardised richness of its FIRST child
(`getchildren(tree, node)[1]`) against the null, which keeps the parent clade's richness
per cell fixed. Positive SOS (the high end) marks cells where the first child is
over-represented; negative SOS (the low end) marks cells where the second child is.
"""
function focuscolors(tree, layout::TreeLayout, node, sos_colormap, contextcolor)
    cmap = Makie.to_colormap(sos_colormap)
    colors = fill(to_color(contextcolor), length(layout))
    for (child, color) in zip(getchildren(tree, node)[1:2], (last(cmap), first(cmap)))
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
[`NodePanel`](@ref) `panel`, and `status`, the text of the label above the tree.
"""
struct NodeExplorer
    figure::Figure
    axis::Axis
    treeplot::TreePlot
    panel::NodePanel
    status::Observable{String}
end

metricvalues(res, metric::Symbol) = getfield(res, metric)
metricvalues(res, metric::AbstractDict) = metric

"""
    nodeexplorer(assemblage, tree, res; kwargs...)

The tree and a [`nodepanel`](@ref) side by side. The tree marks `nodes` coloured by
`metric`; clicking a marker or a branch shows that node in the panel. The two clades
below the node shown are drawn in the end colours of the SOS colour map (see
[`focuscolors`](@ref)), and the rest of the tree in `contextcolor`. `res` is the cached
`NodeMetrics`/`NodeAnalysis`.

Keyword arguments:
- `metric = :gnd`: a field of `res` (`:gnd`, `:rms`, `:ses`, ...) or a Dict of node => value
- `nodes`: the nodes to mark; default every node with an SOS and a finite metric value.
  Pass e.g. `divergent_nodes(res; ...)` to mark only those.
- `node`: the node shown first; default the marked node with the highest metric value
- `treetype = :fan`, `showtips = false`, `colormap = :YlOrRd`, `colorrange = automatic`,
  `markersize = 14`: passed to `treeplot`, as is everything in `treekw = (;)`. The
  branch colours are set by the explorer, so `branchcolor` is overridden.
- `contextcolor = :gray75`: the branches outside the two clades of the node shown
- `panel = (;)`: keyword arguments for `nodepanel!`
- `figure = (;)`: attributes for the `Figure`
- `pickfn = pick`: the picking function (replaced in the tests, where CairoMakie cannot pick)

Returns `(figure, explorer)`; see [`NodeExplorer`](@ref). Clicking needs an interactive
backend (`using GLMakie`); with CairoMakie you get the figure for the initial node.
"""
function nodeexplorer(assemblage, tree, res; metric = :gnd, nodes = automatic,
                      node = automatic, treetype = :fan, showtips = false,
                      colormap = :YlOrRd, colorrange = automatic, markersize = 14,
                      contextcolor = :gray75, treekw = (;), panel = (;), figure = (;),
                      pickfn = pick)
    sos = sosvalues(res)
    vals = metricvalues(res, metric)
    nodes = nodes === automatic ? collect(keys(sos)) : collect(nodes)
    marked = Dict(n => Float64(vals[n]) for n in nodes
                  if haskey(vals, n) && isfinite(vals[n]) && hassos(tree, sos, n))
    isempty(marked) && throw(ArgumentError("None of `nodes` has both an SOS and a finite metric value"))
    node = node === automatic ? argmax(n -> marked[n], keys(marked)) : String(node)
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
    focus(n) = focuscolors(tree, layout, n, sos_colormap, contextcolor)
    tp = treeplot!(ax, tree; treetype, showtips, nodecolor = marked, markersize,
                   colormap, colorrange, merge(treekw, (; branchcolor = focus(node)))...)
    on(n -> (tp.branchcolor = focus(n)), np.node)
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
    return fig, NodeExplorer(fig, ax, tp, np, status)
end
