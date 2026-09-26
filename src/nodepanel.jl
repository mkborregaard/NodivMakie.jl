# Makie version of Nodiv's `plot_node`: a 2x2 grid with the richness of the node's
# clade (top left), its SOS (top right) and the richness of its two descendant clades
# (bottom row). The node is an Observable, so the panel redraws in place when it
# changes; that is what the linked tree-map explorer drives.

"""
    NodePanel

The axes and plots of a [`nodepanel`](@ref). Set `panel.node[] = "name"` to show
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
    axes::Vector{Union{Axis, Nothing}}
    maps::Vector{Any}
    colorbars::Vector{Colorbar}
    sos::AbstractDict
end

sosvalues(res::Nodiv.AbstractNodeResult) = res.sos

"""
    cladecolors(sos_colormap; inset = 0.15) -> (first, second)

The colours standing for a node's two child clades: the high end of the SOS colour map
for the first child (positive SOS: cells where the first child is over-represented) and
the low end for the second, each taken `inset` in from the end.
"""
function cladecolors(sos_colormap; inset = 0.15)
    cmap = Makie.to_colormap(sos_colormap)
    return Makie.interpolated_getindex(cmap, 1.0 - inset),
           Makie.interpolated_getindex(cmap, Float64(inset))
end
sosvalues(d::AbstractDict) = d

"""
    hassos(tree, sos, node)

Whether `node` can be shown in a node panel: it has an SOS vector in `sos` and at least
two children.
"""
hassos(tree, sos, node) = haskey(sos, node) && hasnode(tree, node) &&
                          length(getchildren(tree, node)) >= 2

function occupied(r)
    v = Float64.(r)
    v[v .== 0] .= NaN
    return v
end

# The four maps' values and titles for one node. `richness_of` is Nodiv's
# `clade_richness(assemblage, tree)`, fast enough for clicking through nodes.
function nodepaneldata(richness_of, tree, sos, node)
    ch1, ch2 = [getnodename(tree, c) for c in getchildren(tree, node)[1:2]]
    return (values = [occupied(richness_of(node)), sitevalues(sos[node]),
                      occupied(richness_of(ch1)), occupied(richness_of(ch2))],
            titles = [node, "SOS", ch1, ch2])
end

"""
    nodepanel!(gridposition, assemblage, tree, node, res; kwargs...)

Draw the node panel for `node` (a node name or an `Observable` of one) into a figure
position. `res` is a `NodeAnalysis`/`NodeMetrics` or a Dict of node name => SOS
vector, e.g. a cached result; nothing is recomputed.

Keyword arguments:
- `richness_colormap = Reverse(:Spectral)`, `sos_colormap = :RdYlBu`,
  `sos_colorrange = (-8, 8)`: the colours of `plot_node`
- `colorbars = true`: a colour bar beside each map
- `clademap = true`: draw the map of the node's own clade (top left). With `false` that
  cell of `panel.layout` is left free, as the [`nodeexplorer`](@ref) does for its ordination.
- `titlecolors = true`: colour the child clades' map titles by [`cladecolors`](@ref), as
  the explorer colours their branches (`colorinset` as its `focusinset`)
- `axis = (;)`: attributes for all four axes
- `images = nothing`: species images for the child-clade maps, as a
  [`SpeciesImages`](@ref) or a directory; each gets its widest-ranging species with an
  image in the top-right corner (see [`cladeimages!`](@ref)). Range sizes are taken from
  `assemblage` unless `imageoptions` has a `rangesize`.
- `imageoptions = (;)`: keyword arguments for `cladeimages!`

Returns a [`NodePanel`](@ref).
"""
function nodepanel!(gp, assemblage, tree, node, res;
                    richness_colormap = Reverse(:Spectral), sos_colormap = :RdYlBu,
                    sos_colorrange = (-8, 8), colorbars = true, axis = (;),
                    images = nothing, imageoptions = (;), titlecolors = true,
                    colorinset = 0.15, clademap = true)
    sos = sosvalues(res)
    node = node isa Observable ? node : Observable(String(node))
    hassos(tree, sos, node[]) ||
        throw(ArgumentError("Node $(node[]) has no SOS in the result or fewer than two children"))

    # Only valid nodes reach the maps; an invalid one leaves the panel as it was
    richness_of = clade_richness(assemblage, tree)
    data = Observable(nodepaneldata(richness_of, tree, sos, node[]))
    on(node) do n
        if hassos(tree, sos, n)
            data[] = nodepaneldata(richness_of, tree, sos, n)
        else
            @warn "Node $n has no SOS in the result or fewer than two children; not shown"
        end
    end

    gl = GridLayout(gp)
    locs = sitelocations(assemblage)
    axes, maps, cbs = Union{Axis, Nothing}[], Any[], Colorbar[]
    for (i, (row, col)) in enumerate(((1, 1), (1, 2), (2, 1), (2, 2)))
        if i == 1 && !clademap
            push!(axes, nothing)
            push!(maps, nothing)
            continue
        end
        ax = Axis(gl[row, 2col - 1]; title = lift(d -> d.titles[i], data),
                  autolimitaspect = 1, xgridvisible = false, ygridvisible = false, axis...)
        cm = i == 2 ? (colormap = sos_colormap, colorrange = sos_colorrange) :
             (colormap = richness_colormap,)
        m = sitemap!(ax, lift(d -> d.values[i], data), locs; cm...)
        push!(axes, ax)
        push!(maps, m)
        colorbars && push!(cbs, Colorbar(gl[row, 2col], m; width = 10))
    end
    linkaxes!(filter(!isnothing, axes)...)
    if titlecolors
        axes[3].titlecolor, axes[4].titlecolor = cladecolors(sos_colormap; inset = colorinset)
    end
    np = NodePanel(node, gl, axes, maps, cbs, sos)
    if images !== nothing
        rangesize = get(imageoptions, :rangesize, assemblage)
        cladeimages!(np, tree, images, rangesize;
                     Base.structdiff(imageoptions, NamedTuple{(:rangesize,)})...)
    end
    return np
end

"""
    nodepanel(assemblage, tree, node, res; figure = (;), kwargs...)

[`nodepanel!`](@ref) in a new figure. Returns `(figure, panel)`.
"""
function nodepanel(assemblage, tree, node, res; figure = (;), kwargs...)
    fig = Figure(; size = (900, 800), figure...)
    return fig, nodepanel!(fig[1, 1], assemblage, tree, node, res; kwargs...)
end

"""
    sosmap!(gridposition, assemblage, node, res; kwargs...) -> (axis, map)

The SOS of `node` (a node name or an `Observable` of one) as a map with a colour bar,
from the cached `res` (a `NodeAnalysis`/`NodeMetrics` or a Dict of node name => SOS
vector). The map follows `node`; a node with no SOS in `res` leaves it as it was. With
[`explorertree!`](@ref) this builds explorers with other panels, e.g. the SOS of the
same node in two spaces side by side.

Keyword arguments:
- `title = "SOS"`: the axis title, or a function of the node name giving it
- `colormap = :RdYlBu`, `colorrange = (-8, 8)`: as the SOS map of [`nodepanel!`](@ref)
- `colorbar = true`: a colour bar beside the map
- `axis = (;)`: attributes for the axis
"""
function sosmap!(gp, assemblage, node, res; title = "SOS", colormap = :RdYlBu,
                 colorrange = (-8, 8), colorbar = true, axis = (;))
    sos = sosvalues(res)
    node = node isa Observable ? node : Observable(String(node))
    haskey(sos, node[]) || throw(ArgumentError("Node $(node[]) has no SOS in the result"))
    shown = Observable(node[])
    on(n -> haskey(sos, n) && (shown[] = n), node)
    gl = GridLayout(gp)
    ax = Axis(gl[1, 1]; title = title isa Function ? lift(title, shown) : title,
              autolimitaspect = 1, xgridvisible = false, ygridvisible = false, axis...)
    m = sitemap!(ax, lift(n -> sitevalues(sos[n]), shown), sitelocations(assemblage);
                 colormap, colorrange)
    colorbar && Colorbar(gl[1, 2], m; width = 10)
    return ax, m
end
