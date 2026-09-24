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
- `axes`: the four map axes (clade, SOS, child 1, child 2), with linked limits
- `maps`: the four `SiteMap` plots, in the same order
- `colorbars`: their colour bars (empty if `colorbars = false`)
"""
struct NodePanel
    node::Observable{String}
    layout::GridLayout
    axes::Vector{Axis}
    maps::Vector{Any}
    colorbars::Vector{Colorbar}
end

sosvalues(res::Union{Nodiv.NodeAnalysis, Nodiv.NodeMetrics}) = res.sos
sosvalues(d::AbstractDict) = d

"""
    hassos(tree, sos, node)

Whether `node` can be shown in a node panel: it has an SOS vector in `sos` and at least
two children.
"""
hassos(tree, sos, node) = haskey(sos, node) && hasnode(tree, node) &&
                          length(getchildren(tree, node)) >= 2

# Clade richness for interactive use. `richness(get_clade(...))` on a species view of a
# large assemblage takes ~0.4 s (birds, 18k cells), too slow for clicking. Here the
# presence matrix is stored once as sites x species, so a clade's richness is a sum over
# its species' columns. Same result as `richness(get_clade(assemblage, tree, node))`.
struct CladeRichness
    presence::SparseMatrixCSC{Bool, Int}   # sites x species
    index::Dict{String, Int}               # species name => column
end

function CladeRichness(assemblage)
    occ = sparse(EcoBase.occurrences(assemblage) .> 0)            # species x sites
    names = String.(EcoBase.thingnames(assemblage))
    return CladeRichness(permutedims(occ), Dict(n => i for (i, n) in enumerate(names)))
end

function (cr::CladeRichness)(tree, node)
    cols = [cr.index[sp] for sp in nodespecies(tree, node) if haskey(cr.index, sp)]
    r = zeros(Int, size(cr.presence, 1))
    rows, vals = rowvals(cr.presence), nonzeros(cr.presence)
    for j in cols, k in nzrange(cr.presence, j)
        vals[k] && (r[rows[k]] += 1)
    end
    return r
end

function occupied(r)
    v = Float64.(r)
    v[v .== 0] .= NaN
    return v
end

# The four maps' values and titles for one node
function nodepaneldata(cr::CladeRichness, tree, sos, node)
    ch1, ch2 = [getnodename(tree, c) for c in getchildren(tree, node)[1:2]]
    return (values = [occupied(cr(tree, node)), sitevalues(sos[node]),
                      occupied(cr(tree, ch1)), occupied(cr(tree, ch2))],
            titles = [node, "SOS", ch1, ch2])
end

"""
    nodepanel!(gridposition, assemblage, tree, node, res; kwargs...)

Draw the node panel for `node` (a node name or an `Observable` of one) into a figure
position. `res` is a `NodeAnalysis`/`NodeMetrics` or a Dict of node name => SOS
vector, e.g. the cached `res_e`/`res_g`; nothing is recomputed.

Keyword arguments:
- `richness_colormap = Reverse(:Spectral)`, `sos_colormap = :RdYlBu`,
  `sos_colorrange = (-8, 8)`: the colours of `plot_node`
- `colorbars = true`: a colour bar beside each map
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
                    images = nothing, imageoptions = (;))
    sos = sosvalues(res)
    node = node isa Observable ? node : Observable(String(node))
    hassos(tree, sos, node[]) ||
        throw(ArgumentError("Node $(node[]) has no SOS in the result or fewer than two children"))

    # Only valid nodes reach the maps; an invalid one leaves the panel as it was
    cr = CladeRichness(assemblage)
    data = Observable(nodepaneldata(cr, tree, sos, node[]))
    on(node) do n
        if hassos(tree, sos, n)
            data[] = nodepaneldata(cr, tree, sos, n)
        else
            @warn "Node $n has no SOS in the result or fewer than two children; not shown"
        end
    end

    gl = GridLayout(gp)
    locs = sitelocations(assemblage)
    axes, maps, cbs = Axis[], Any[], Colorbar[]
    for (i, (row, col)) in enumerate(((1, 1), (1, 2), (2, 1), (2, 2)))
        ax = Axis(gl[row, 2col - 1]; title = lift(d -> d.titles[i], data),
                  autolimitaspect = 1, xgridvisible = false, ygridvisible = false, axis...)
        cm = i == 2 ? (colormap = sos_colormap, colorrange = sos_colorrange) :
             (colormap = richness_colormap,)
        m = sitemap!(ax, lift(d -> d.values[i], data), locs; cm...)
        push!(axes, ax)
        push!(maps, m)
        colorbars && push!(cbs, Colorbar(gl[row, 2col], m; width = 10))
    end
    linkaxes!(axes...)
    np = NodePanel(node, gl, axes, maps, cbs)
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
