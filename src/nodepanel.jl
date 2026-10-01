# Makie version of Nodiv's `plot_node`: a 2x2 grid with the richness of the node's
# clade (top left), its SOS (top right) and the richness of its two descendant clades
# (bottom row). The node is an Observable, so the panel redraws in place when it
# changes; that is what the linked tree-map explorer drives.

"""
    NodePanel

The axes and plots of a [`node_panel`](@ref). Set `panel.node[] = "name"` to show
another node.

- `node`: `Observable{String}` with the node shown
- `layout`: the `GridLayout` holding the panel
- `axes`: the four map axes (clade, SOS, child 1, child 2), with linked limits. With
  `clademap = false` the first is `nothing`, so the others keep their positions.
- `maps`: the four `SiteMap` plots, in the same order
- `colorbars`: their colour bars (empty if `colorbars = false`)
- `sos`: the SOS vectors the panel shows, a Dict of node name => SOS vector
"""
struct NodePanel
    node::Observable{String}
    layout::GridLayout
    axes::Vector{Union{Axis,Nothing}}
    maps::Vector{Any}
    colorbars::Vector{Colorbar}
    sos::AbstractDict
end

# The SOS vectors of a result, or a Dict of node name => SOS vector as it is
_sos_values(res::Nodiv.AbstractNodeResult) = res.sos
_sos_values(d::AbstractDict) = d

# A node given as a name or as an `Observable` of one
_node_observable(node::Observable) = node
_node_observable(node) = Observable(String(node))

"""
    clade_colors(sos_colormap; inset=0.15) -> (first, second)

The colours standing for a node's two child clades: the high end of the SOS colour map
for the first child (positive SOS: cells where the first child is over-represented) and
the low end for the second, each taken `inset` in from the end.
"""
function clade_colors(sos_colormap; inset=0.15)
    return _colormap_at(sos_colormap, 1 - inset), _colormap_at(sos_colormap, inset)
end

"""
    has_sos(tree, sos, node)

Whether `node` can be shown in a node panel: it has an SOS vector in `sos` and at least
two children.
"""
function has_sos(tree, sos, node)
    return haskey(sos, node) && hasnode(tree, node) && length(getchildren(tree, node)) >= 2
end

# The four maps' values and titles for one node. `richness_of` is Nodiv's
# `clade_richness(assemblage, tree)`, fast enough for clicking through nodes.
function _node_panel_data(richness_of, tree, sos, node)
    ch1, ch2 = [getnodename(tree, c) for c in getchildren(tree, node)[1:2]]
    return (
        values=[
            _occupied(richness_of(node)),
            _site_values(sos[node]),
            _occupied(richness_of(ch1)),
            _occupied(richness_of(ch2)),
        ],
        titles=[node, "SOS", ch1, ch2],
    )
end

"""
    node_panel!(gridposition, assemblage, tree, node, res; kwargs...)

Draw the node panel for `node` (a node name or an `Observable` of one) into a figure
position. `res` is a `NodeAnalysis`/`NodeMetrics` or a Dict of node name => SOS
vector, e.g. a cached result; nothing is recomputed.

Keyword arguments:
- `richness_colormap = Reverse(:Spectral)`, `sos_colormap = :RdYlBu`,
  `sos_colorrange = (-8, 8)`: the colours of `plot_node`
- `sos_empty_color = :lightgray`: the sites with no SOS on the SOS map (see
  [`sitemap`](@ref)'s `empty_color`)
- `colorbars = true`: a colour bar beside each map
- `clademap = true`: draw the map of the node's own clade (top left). With `false` that
  cell of `panel.layout` is left free, as the [`node_explorer`](@ref) does for its
  ordination.
- `titlecolors = true`: colour the child clades' map titles by [`clade_colors`](@ref), as
  the explorer colours their branches (`colorinset` as its `focusinset`)
- `axis = (;)`: attributes for all four axes
- `images = nothing`: species images for the child-clade maps, as a
  [`SpeciesImages`](@ref) or a directory; each gets its widest-ranging species with an
  image in the top-right corner (see [`clade_images!`](@ref)). Range sizes are taken from
  `assemblage` unless `imageoptions` has a `rangesize`.
- `imageoptions = (;)`: keyword arguments for `clade_images!`

Returns a [`NodePanel`](@ref).
"""
function node_panel!(
    gp,
    assemblage,
    tree,
    node,
    res;
    richness_colormap=Reverse(:Spectral),
    sos_colormap=:RdYlBu,
    sos_colorrange=(-8, 8),
    sos_empty_color=:lightgray,
    colorbars=true,
    axis=(;),
    images=nothing,
    imageoptions=(;),
    titlecolors=true,
    colorinset=0.15,
    clademap=true,
)
    sos = _sos_values(res)
    node = _node_observable(node)
    has_sos(tree, sos, node[]) || throw(
        ArgumentError("Node $(node[]) has no SOS in the result or fewer than two children"),
    )

    # Only valid nodes reach the maps; an invalid one leaves the panel as it was
    richness_of = clade_richness(assemblage, tree)
    data = Observable(_node_panel_data(richness_of, tree, sos, node[]))
    on(node) do n
        if has_sos(tree, sos, n)
            data[] = _node_panel_data(richness_of, tree, sos, n)
        else
            @warn "Node $n has no SOS in the result or fewer than two children; not shown"
        end
    end

    gl = GridLayout(gp)
    locs = _site_locations(assemblage)
    axes, maps, cbs = Union{Axis,Nothing}[], Any[], Colorbar[]
    for (i, (row, col)) in enumerate(((1, 1), (1, 2), (2, 1), (2, 2)))
        if i == 1 && !clademap
            push!(axes, nothing)
            push!(maps, nothing)
            continue
        end
        ax = Axis(
            gl[row, 2col - 1];
            title=lift(d -> d.titles[i], data),
            _EQUAL_SCALES...,
            xgridvisible=false,
            ygridvisible=false,
            axis...,
        )
        if i == 2
            cm = (
                colormap=sos_colormap,
                colorrange=sos_colorrange,
                empty_color=sos_empty_color,
            )
        else
            cm = (colormap=richness_colormap,)
        end
        m = sitemap!(ax, lift(d -> d.values[i], data), locs; cm...)
        push!(axes, ax)
        push!(maps, m)
        colorbars && push!(cbs, Colorbar(gl[row, 2col], m; width=10))
    end
    linkaxes!(filter(!isnothing, axes)...)
    if titlecolors
        axes[3].titlecolor, axes[4].titlecolor = clade_colors(
            sos_colormap; inset=colorinset
        )
    end
    np = NodePanel(node, gl, axes, maps, cbs, sos)
    if images !== nothing
        rangesize, opts = _image_options(imageoptions, assemblage)
        clade_images!(np, tree, images, rangesize; opts...)
    end
    return np
end

"""
    node_panel(assemblage, tree, node, res; figure=(;), kwargs...)

[`node_panel!`](@ref) in a new figure. Returns `(figure, panel)`.
"""
function node_panel(assemblage, tree, node, res; figure=(;), kwargs...)
    fig = Figure(; size=(900, 800), figure...)
    return fig, node_panel!(fig[1, 1], assemblage, tree, node, res; kwargs...)
end

"""
    sos_map!(gridposition, assemblage, node, res; kwargs...) -> (axis, map)

The SOS of `node` (a node name or an `Observable` of one) as a map with a colour bar,
from the cached `res` (a `NodeAnalysis`/`NodeMetrics` or a Dict of node name => SOS
vector). The map follows `node`; a node with no SOS in `res` leaves it as it was. With
[`explorer_tree!`](@ref) this builds explorers with other panels, e.g. the SOS of the
same node in two spaces side by side.

Keyword arguments:
- `title = "SOS"`: the axis title, or a function of the node name giving it
- `colormap = :RdYlBu`, `colorrange = (-8, 8)`, `empty_color = :lightgray`: as the SOS
  map of [`node_panel!`](@ref)
- `colorbar = true`: a colour bar beside the map
- `axis = (;)`: attributes for the axis
"""
function sos_map!(
    gp,
    assemblage,
    node,
    res;
    title="SOS",
    colormap=:RdYlBu,
    colorrange=(-8, 8),
    empty_color=:lightgray,
    colorbar=true,
    axis=(;),
)
    sos = _sos_values(res)
    node = _node_observable(node)
    haskey(sos, node[]) || throw(ArgumentError("Node $(node[]) has no SOS in the result"))
    shown = Observable(node[])
    on(n -> haskey(sos, n) && (shown[] = n), node)
    gl = GridLayout(gp)
    ax = Axis(
        gl[1, 1];
        title=title isa Function ? lift(title, shown) : title,
        _EQUAL_SCALES...,
        xgridvisible=false,
        ygridvisible=false,
        axis...,
    )
    m = sitemap!(
        ax,
        lift(n -> _site_values(sos[n]), shown),
        _site_locations(assemblage);
        colormap,
        colorrange,
        empty_color,
    )
    colorbar && Colorbar(gl[1, 2], m; width=10)
    return ax, m
end
