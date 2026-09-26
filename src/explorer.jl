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
    onnodeclick(f, ax, plot; range = 10)

Call `f(nodename)` when a node of `plot` in `ax` is left-clicked: a node marker or branch
of a `treeplot`, or a point of an [`ordinationplot`](@ref). `range` is the pick radius in
pixels. Needs an interactive backend (GLMakie, WGLMakie).
"""
onnodeclick(f, ax::Axis, p::Union{TreePlot, OrdinationPlot}; range = 10) =
    _onnodeclick(f, ax, p, pick; range)

# `pickfn` is the backend's `pick`; the tests swap in a stand-in, as CairoMakie cannot pick
function _onnodeclick(f, ax::Axis, p, pickfn; range = 10)
    scene = ax.scene
    return on(events(scene).mousebutton; priority = 2) do ev
        (ev.button == Mouse.left && ev.action == Mouse.press) || return Consume(false)
        Makie.is_mouseinside(scene) || return Consume(false)
        n = nodeat(p, pickfn(scene, events(scene).mouseposition[], range)...)
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
    high, low = cladecolors(sos_colormap; inset)
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
    ExplorerTree

The tree side of an explorer, from [`explorertree!`](@ref): `node`, the `Observable` with
the node shown, `layout`, the `GridLayout` holding the label, tree and colour bar, the
tree `axis` and `treeplot`, `status`, the text of the label above the tree, and `images`,
the species images around the tree ([`TreeImages`](@ref)) or `nothing`.
"""
struct ExplorerTree
    node::Observable{String}
    layout::GridLayout
    axis::Axis
    treeplot::TreePlot
    status::Observable{String}
    images::Union{TreeImages, Nothing}
end

"""
    explorertree!(gridposition, tree, node, marked; kwargs...)

The tree side of an explorer, for building explorers with other panels. Draws the tree
with markers on the nodes of `marked` (a Dict of node name => value), coloured by value,
a label above it with the node shown and a colour bar below. `node` is an
`Observable{String}` with the node shown: the two clades below it are drawn in the SOS
colours (see [`focuscolors`](@ref)), and clicking a node marker or a branch sets it. Any
panel that follows `node` then updates with the tree; see [`sosmap!`](@ref) for one.

Keyword arguments:
- `selectable = n -> true`: whether a clicked node can be shown. A node that cannot is
  reported in the label as `"\$n: \$unselectable"` and not shown.
- `unselectable = "no SOS (a tip, or not analysed)"`
- `values = marked`: node name => value, shown in the label and the hover labels; pass
  every node's value to label unmarked nodes too
- `label = "value"`: the name of the value, in the label and on the colour bar
- `focuscolormap = :RdYlBu`: the SOS colour map the two clades take their colours from;
  `focusinset` and `contextcolor` as for [`nodeexplorer`](@ref)
- `treetype`, `showtips`, `colormap`, `colorrange`, `markersize`, `strokewidth`,
  `strokecolor`, `treekw`, `pickfn`: as for [`nodeexplorer`](@ref)
- `images = nothing`: species images around the tree, as a [`SpeciesImages`](@ref) or a
  directory path, with range sizes from `rangesize` (e.g. the assemblage) and the keyword
  arguments `imageoptions = (;)` for [`treeimages!`](@ref)

The hover labels (the node, its number of species and its value) show once the figure
has a `DataInspector`. Returns an [`ExplorerTree`](@ref).
"""
function explorertree!(gp, tree, node::Observable{String}, marked::AbstractDict;
                       selectable = n -> true,
                       unselectable = "no SOS (a tip, or not analysed)",
                       values = marked, label = "value", treetype = :fan,
                       showtips = false, colormap = :YlOrRd, colorrange = automatic,
                       markersize = 14, strokewidth = 1, strokecolor = :gray20,
                       contextcolor = :gray75, focuscolormap = :RdYlBu, focusinset = 0.15,
                       images = nothing, imageoptions = (;), rangesize = nothing,
                       treekw = (;), pickfn = pick)
    gl = GridLayout(gp)
    status = Observable("")
    Label(gl[1, 1], status; tellwidth = false, halign = :left)
    ax = Axis(gl[2, 1]; autolimitaspect = treetype === :fan ? 1 : nothing)
    hidedecorations!(ax)
    hidespines!(ax)

    # the node shown: its two clades in the SOS colours, the rest of the tree greyed
    layout = treelayout(tree)
    focus(n) = focuscolors(tree, layout, n, focuscolormap, contextcolor; inset = focusinset)
    # markers only on the marked nodes: a NaN colour hides a marker's fill, not its outline
    tp = treeplot!(ax, tree; treetype, showtips, nodecolor = marked,
                   shownodes = collect(keys(marked)), markersize,
                   strokewidth, strokecolor, colormap, colorrange,
                   merge(treekw, (; branchcolor = focus(node[])))...)
    on(n -> (tp.branchcolor = focus(n)), node)
    ti = nothing
    if images !== nothing
        rangesize === nothing &&
            throw(ArgumentError("Species images around the tree need `rangesize`, e.g. the assemblage"))
        images isa AbstractString && (images = SpeciesImages(images))
        ti = treeimages!(ax, tp, images, rangesize; imageoptions...)
    end
    Colorbar(gl[3, 1], tp; vertical = false, flipaxis = false, label,
             tellheight = true, width = Relative(0.6))

    describe(n) = haskey(values, n) && values[n] isa Real ? "$n   $label = $(round(values[n]; sigdigits = 3))" : n
    status[] = describe(node[])
    on(n -> (status[] = describe(n)), node)
    _onnodeclick(ax, tp, pickfn) do n
        if selectable(n)
            node[] = n
        else
            status[] = "$n: $unselectable"
        end
    end
    # hover labels: the node, its size, and its value
    function hover(n)
        k = tp.tree_layout[].index[n]
        lab = tp.tree_layout[].isleaf[k] ? n : "$n  ($(tp.clade_sizes[][k]) species)"
        haskey(values, n) && values[n] isa Real && isfinite(values[n]) &&
            (lab *= "\n$label = $(round(values[n]; sigdigits = 3))")
        return lab
    end
    tp.hoverlabel = hover
    return ExplorerTree(node, gl, ax, tp, status, ti)
end

"""
    NodeExplorer

The parts of a [`nodeexplorer`](@ref): `figure`, the tree `axis` and `treeplot`, the
[`NodePanel`](@ref) `panel`, `status`, the text of the label above the tree,
`images`, the species images around the tree ([`TreeImages`](@ref)) or `nothing`,
`inspector`, the `DataInspector` showing hover labels (or `nothing`), and
`ordination` and `ordinationaxis`, the [`ordinationplot`](@ref) of the marked nodes and
its axis (both `nothing` when the panel shows the clade map instead).
"""
struct NodeExplorer
    figure::Figure
    axis::Axis
    treeplot::TreePlot
    panel::NodePanel
    status::Observable{String}
    images::Union{TreeImages, Nothing}
    inspector::Union{DataInspector, Nothing}
    ordination::Union{OrdinationPlot, Nothing}
    ordinationaxis::Union{Axis, Nothing}
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

In the panel, the top-left cell shows an ordination of the marked nodes by the
similarity of their SOS maps ([`sosordination`](@ref)) instead of the map of the node's
own clade. Its points are coloured like the tree's markers. Clicking a point shows that
node, just like clicking it on the tree, and the node shown has a ring. With fewer than
three marked nodes, or `ordination = false`, the cell shows the clade map instead.

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
- `images = nothing`: species images, as a [`SpeciesImages`](@ref) or a directory path:
  drawn around the tree (see [`treeimages!`](@ref)) and in the corner of the two
  child-clade maps at the same size (see [`cladeimages!`](@ref)). Range sizes are taken
  from `assemblage` unless `imageoptions` has a `rangesize`.
- `imageoptions = (;)`: keyword arguments for `treeimages!`; `shape`, `fit`,
  `whitebackground` and `clip` apply to the map images too
- `ordination = true`: the ordination of the marked nodes in place of the clade map
- `ordinationkw = (;)`: keyword arguments for [`sosordination`](@ref), e.g. the
  `minoverlap` of `sos_distances` (default 3), which should be set for each space
- `panel = (;)`: keyword arguments for `nodepanel!`
- `figure = (;)`: attributes for the `Figure`
- `inspector = true`: hovering shows a label (a `DataInspector`): on the tree and the
  ordination the node (for a branch, the node it leads to) with its number of species and
  `metric` value, on the maps the cell's value. Needs an interactive backend.
- `pickfn = pick`: the picking function (replaced in the tests, where CairoMakie cannot pick)

Returns `(figure, explorer)`; see [`NodeExplorer`](@ref). Clicking needs an interactive
backend (`using GLMakie`); with CairoMakie you get the figure for the initial node.
The tree side is [`explorertree!`](@ref), which builds explorers with other panels.
"""
function nodeexplorer(assemblage, tree, res; metric = defaultmetric(res),
                      nodes = automatic, node = automatic, treetype = :fan,
                      showtips = false, colormap = :YlOrRd, colorrange = automatic,
                      markersize = 14, strokewidth = 1, strokecolor = :gray20,
                      contextcolor = :gray75, focusinset = 0.15, images = nothing,
                      imageoptions = (;), treekw = (;), ordination = true,
                      ordinationkw = (;), panel = (;), figure = (;),
                      inspector = true, pickfn = pick)
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
    showordination = ordination && length(marked) >= 3
    np = nodepanel!(fig[1, 2], assemblage, tree, node, res;
                    merge((; colorinset = focusinset, clademap = !showordination), panel)...)
    rangesize = get(imageoptions, :rangesize, assemblage)
    opts = Base.structdiff(imageoptions, NamedTuple{(:rangesize,)})
    images isa AbstractString && (images = SpeciesImages(images))
    et = explorertree!(fig[1, 1], tree, np.node, marked; values = vals, label,
                       selectable = n -> hassos(tree, sos, n), treetype, showtips,
                       colormap, colorrange, markersize, strokewidth, strokecolor,
                       contextcolor, focuscolormap = np.maps[2].colormap[], focusinset,
                       images, imageoptions = opts, rangesize, treekw, pickfn)
    colsize!(fig.layout, 1, Relative(0.45))
    tp = et.treeplot
    if images !== nothing
        shared = NamedTuple(k => v for (k, v) in pairs(opts) if k in (:shape, :fit, :whitebackground, :clip))
        cladeimages!(np, tree, images, rangesize; pixelsize = et.images.pixelsize, shared...)
    end

    # the marked nodes by SOS similarity, coloured like the tree's markers, with the tree's
    # hover labels; a click on a point selects its node, and the node shown has a ring
    op, oax = nothing, nothing
    if showordination
        ord = sosordination(res, sort!(collect(keys(marked))); ordinationkw...)
        oax = Axis(np.layout[1, 1]; title = "SOS similarity ($(length(ord.nodes)) nodes)",
                   xlabel = "MDS axis 1", ylabel = "MDS axis 2", autolimitaspect = 1,
                   xgridvisible = false, ygridvisible = false)
        op = ordinationplot!(oax, ord; nodecolor = marked, colormap,
                             colorrange = tp.joint_colorrange[], markersize, strokewidth,
                             strokecolor, selected = node, hoverlabel = tp.hoverlabel[])
        on(n -> (op.selected = n), np.node)
        _onnodeclick(n -> (np.node[] = n), oax, op, pickfn)
    end
    di = inspector ? DataInspector(fig) : nothing
    return fig, NodeExplorer(fig, et.axis, tp, np, et.status, et.images, di, op, oax)
end
