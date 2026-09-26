# Plots of Nodiv's `sos_ordination`, the ordination of nodes by the similarity of their SOS
# maps: its points, picked and clicked like tree nodes, and its eigenvalues.

"""
    ordinationplot(ordination; kwargs...)

The nodes of an `SOSOrdination` (from Nodiv's `sos_ordination`) as points on its first two
axes. Like `treeplot`'s, the per-node values (`nodecolor`, `nodelabels`) are a Dict of node name =>
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

The eigenvalues of the axes of an `SOSOrdination` (from Nodiv's `sos_ordination`) as bars,
largest first. They tell whether the first two axes, the ones [`ordinationplot`](@ref)
shows, carry the structure: if axes 3 and up are about as large, the 2-D picture is a
projection. Fit the ordination with more than two axes for this, e.g.
`sos_ordination(D, nodes; maxoutdim = 10)`.
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
