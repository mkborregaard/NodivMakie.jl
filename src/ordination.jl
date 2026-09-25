# Ordination of nodes by the similarity of their SOS maps: classical MDS on Nodiv's
# `sos_distances`, and a plot of it whose points are picked and clicked like tree nodes.

"""
    SOSOrdination

The result of [`sosordination`](@ref):
- `nodes`: the node names, in the order of the columns of `coords`
- `coords`: the MDS coordinates, one column per node and one row per axis
- `eigenvalues`: the eigenvalues of the axes, largest first
- `distances`: the SOS distance matrix the ordination was fitted to
"""
struct SOSOrdination
    nodes::Vector{String}
    coords::Matrix{Float64}
    eigenvalues::Vector{Float64}
    distances::Matrix{Float64}
end

sosdistances(res::Union{Nodiv.NodeAnalysis, Nodiv.NodeMetrics}, nodes; kwargs...) =
    sos_distances(res, nodes; kwargs...)
sosdistances(sos::AbstractDict, nodes; kwargs...) = sos_distances([sos[n] for n in nodes]; kwargs...)

"""
    sosordination(res, nodes; maxoutdim = 2, kwargs...)
    sosordination(D::AbstractMatrix, nodes; maxoutdim = 2)

Classical (metric) MDS of `nodes` by the similarity of their SOS maps. The distances are
Nodiv's `sos_distances(res, nodes; kwargs...)`: 1 - |r| between two nodes' SOS over the
cells both occupy, with pairs sharing fewer than `minoverlap` cells (default 3) put at
distance 1. Set `minoverlap` for each space: a floor of about 8 cells suits a
geographic grid of thousands of cells, and a smaller one an environmental space of tens of
bins. `res` is a `NodeAnalysis`/`NodeMetrics` or a Dict of node name => SOS vector;
nothing is recomputed. Or pass the distance matrix `D` of `nodes` itself, e.g. one
computed once with `sos_distances` and shared with other analyses of the same nodes.

`maxoutdim` is the number of axes (at most `length(nodes) - 1`). The `eigenvalues` of an
ordination with more axes (say 10) tell whether two show the structure. Where the nodes
are mostly unrelated in SOS pattern, the distances are all near 1 and the points form a
ring. That is the finding, not a failure of the method.

Returns an [`SOSOrdination`](@ref); plot it with [`ordinationplot`](@ref), and its
eigenvalues with [`eigenvalueplot`](@ref).
"""
function sosordination(res, nodes; maxoutdim = 2, kwargs...)
    nodes = checkordnodes(nodes)
    return sosordination(sosdistances(res, nodes; kwargs...), nodes; maxoutdim)
end

function sosordination(D::AbstractMatrix{<:Real}, nodes; maxoutdim = 2)
    nodes = checkordnodes(nodes)
    size(D) == (length(nodes), length(nodes)) ||
        throw(ArgumentError("The distance matrix must be $(length(nodes)) x $(length(nodes)), one row and column per node; got $(join(size(D), " x "))"))
    D = Matrix{Float64}(D)
    M = MultivariateStats.fit(MultivariateStats.MDS, D; distances = true, maxoutdim)
    return SOSOrdination(nodes, MultivariateStats.predict(M), MultivariateStats.eigvals(M), D)
end

function checkordnodes(nodes)
    nodes = String[n for n in nodes]
    length(nodes) >= 3 ||
        throw(ArgumentError("An ordination needs at least 3 nodes; got $(length(nodes))"))
    return nodes
end

"""
    ordinationplot(ordination; kwargs...)

The nodes of an [`SOSOrdination`](@ref) as points on its first two axes. Like
`treeplot`'s, the per-node values (`nodecolor`, `nodelabels`) are a Dict of node name =>
value (nodes missing from it get `nan_color`, transparent by default), a vector in
`ordination.nodes` order, or a function of the node name. `selected` marks one node with
a ring. [`nodeat`](@ref) and [`onnodeclick`](@ref) work on it as on a tree, which is how
the [`nodeexplorer`](@ref) selects a node from its ordination.
"""
@recipe OrdinationPlot (ordination,) begin
    "Point colour: a colour, or per-node values mapped through `colormap`."
    nodecolor = :gray40
    "Point size."
    markersize = 10
    marker = :circle
    strokewidth = 1
    strokecolor = :gray20
    "Text at the points: `true` for the node names, or per-node strings."
    nodelabels = nothing
    nodelabelsize = 8
    nodelabelcolor = :black
    "The node marked with a ring, or `nothing`. A node that is not in the ordination gets no ring."
    selected = nothing
    "The colour of the ring."
    selectedcolor = :black
    """
    The text shown on hovering over a point (with a `DataInspector`): a function of the node
    name. `automatic` gives the name.
    """
    hoverlabel = automatic
    Makie.mixin_colormap_attributes()...
    Makie.mixin_generic_plot_attributes()...
end

Makie.convert_arguments(::Type{<:OrdinationPlot}, o::SOSOrdination) = (o,)

# One value per node of the ordination, `missing` where the spec has none
ordvalues(x::AbstractDict, nodes) = Any[get(x, n, missing) for n in nodes]
ordvalues(f::Function, nodes) = Any[f(n) for n in nodes]
function ordvalues(x::AbstractVector, nodes)
    length(x) == length(nodes) ||
        throw(ArgumentError("A per-node vector must have one value per node ($(length(nodes))); got $(length(x))"))
    return Vector{Any}(x)
end

function Makie.plot!(p::OrdinationPlot)
    map!(p, :ordination, :points) do o
        y = size(o.coords, 1) >= 2 ? o.coords[2, :] : zeros(length(o.nodes))
        return Point2d.(o.coords[1, :], y)
    end
    map!(p, [:ordination, :nodecolor, :nan_color], :point_colors) do o, nc, nanc
        return iscolorspec(nc) ? nc : colorvalues(ordvalues(nc, o.nodes), nanc)
    end
    map!(p, [:ordination, :nodelabels], :label_texts) do o, lab
        (lab === nothing || lab === false) && return fill("", length(o.nodes))
        lab === true && return copy(o.nodes)
        return String[v isa Union{Missing, Nothing} ? "" : string(v) for v in ordvalues(lab, o.nodes)]
    end
    map!(p, [:ordination, :points, :selected], :selected_points) do o, pts, sel
        i = sel === nothing ? nothing : findfirst(==(sel), o.nodes)
        return i === nothing ? Point2d[] : [pts[i]]
    end
    map!(ms -> 1.8 * ms, p, :markersize, :ring_size)

    function hovertext(plt, i, _)
        n = nodeat(p, plt, i)
        n === nothing && return ""
        f = p.hoverlabel[]
        return f === automatic ? n : string(f(n))
    end

    scatter!(p, p.points; color = p.point_colors, markersize = p.markersize,
             marker = p.marker, strokewidth = p.strokewidth, strokecolor = p.strokecolor,
             colormap = p.colormap, colorscale = p.colorscale, colorrange = p.colorrange,
             lowclip = p.lowclip, highclip = p.highclip, nan_color = p.nan_color,
             alpha = p.alpha, inspector_label = hovertext)
    scatter!(p, p.selected_points; color = :transparent, markersize = p.ring_size,
             strokewidth = 2, strokecolor = p.selectedcolor, inspectable = false)
    text!(p, p.points; text = p.label_texts, fontsize = p.nodelabelsize,
          color = p.nodelabelcolor, align = (:center, :bottom), offset = (0, 4),
          inspectable = false)
    return p
end

# Axis defaults when `ordinationplot` creates the axis: MDS distances are only faithful
# with equal scales on both axes
Makie.preferred_axis_attributes(::Type{Axis}, ::OrdinationPlot) =
    (; autolimitaspect = 1, xlabel = "MDS axis 1", ylabel = "MDS axis 2",
     xgridvisible = false, ygridvisible = false)

# `Colorbar(fig[1, 2], p)`: the colours of the points
Makie.extract_colormap(p::OrdinationPlot) = Makie.extract_colormap(p.plots[1])

"""
    eigenvalueplot(ordination; kwargs...)

The eigenvalues of the axes of an [`SOSOrdination`](@ref) as bars, largest first. They
tell whether the first two axes, the ones [`ordinationplot`](@ref) shows, carry the
structure: if axes 3 and up are about as large, the 2-D picture is a projection. Fit the
ordination with more than two axes for this, e.g. `sosordination(D, nodes; maxoutdim = 10)`.
"""
@recipe EigenvaluePlot (ordination,) begin
    "Bar colour."
    color = :gray40
    Makie.mixin_generic_plot_attributes()...
end

Makie.convert_arguments(::Type{<:EigenvaluePlot}, o::SOSOrdination) = (o,)

function Makie.plot!(p::EigenvaluePlot)
    map!(o -> Point2d.(eachindex(o.eigenvalues), o.eigenvalues), p, :ordination, :bars)
    barplot!(p, p.bars; color = p.color)
    return p
end

Makie.preferred_axis_attributes(::Type{Axis}, ::EigenvaluePlot) =
    (; xlabel = "MDS axis", ylabel = "eigenvalue", xgridvisible = false)

"""
    nodeat(ordinationplot, plot, index)

The node name under a pick result `(plot, index)` on an [`ordinationplot`](@ref): a
point gives its node. Anything else gives `nothing`.
"""
function nodeat(op::OrdinationPlot, plt, idx)
    (plt === nothing || idx === nothing || idx < 1) && return nothing
    plt === op.plots[1] || return nothing          # the points, not the ring or labels
    nodes = op.ordination[].nodes
    return idx <= length(nodes) ? nodes[idx] : nothing
end
