# What NodivMakie uses from Makie, GridLayoutBase and Phylo that they do not mark as
# public. It is all here, so a change upstream needs fixing in one place; the tests cover
# each use. Everything else comes from their public API.

# Phylo's node positions: heights, depths and the nodes in plotting order, as its Plots
# recipe draws them. `_findxy` is internal to Phylo.
_phylo_xy(tree) = Phylo._findxy(tree)

# The colour at `x` (0 to 1) along a colour map, interpolated between its colours.
# `interpolated_getindex` is how Makie samples its colour maps, but it is not public.
_colormap_at(colormap, x) = Makie.interpolated_getindex(to_colormap(colormap), Float64(x))

# A hierarchical clustering as the nodes of Makie's `dendrogram`, with the merges at their
# heights. Makie's `dendrogram` has no public method taking a `Hclust` or merge heights;
# `hcl_nodes` converts one, but it is not public.
_dendrogram_nodes(hc) = Makie.hcl_nodes(hc; useheight=true)

# Show `text` in the `DataInspector`'s tooltip, at the cursor, from a plot's
# `inspector_hover`. This is how Makie's documentation writes a custom tooltip
# ("Extending the DataInspector"), but `update_tooltip_alignment!` is not public.
function _show_tooltip!(inspector, text)
    pos = Point2f(mouseposition_px(inspector.root))
    Makie.update_tooltip_alignment!(inspector, pos; text)
    return nothing
end

# The layout an axis is placed in, and the rows and columns it spans there, for putting a
# column beside it. GridLayoutBase has no public way to ask a block where it is.
function _grid_position(ax)
    gc = GridLayoutBase.gridcontent(ax)
    return gc.parent, gc.span.rows, gc.span.cols
end

# `Colorbar(fig[1, 2], p)` for the recipes. `extract_colormap` is how Makie's own recipes
# tell a colour bar which colours they use, but it is undocumented.

# The node colours if they are numeric, else the branches'
function Makie.extract_colormap(p::TreePlot)
    for child in Iterators.reverse(p.plots)
        child isa Union{Scatter,Lines} || continue
        cm = Makie.extract_colormap(child)
        cm isa Makie.ColorMapping && return cm
    end
    return nothing
end

# The colours of the points
Makie.extract_colormap(p::OrdinationPlot) = Makie.extract_colormap(p.plots[1])

# `Legend` and `axislegend` entries for a tree are its node groups (if any), not the tree
# itself. `get_plots` is how Makie's legends find a plot's entries, but it is undocumented.
function Makie.get_plots(p::TreePlot)
    return filter(c -> c isa Scatter && haskey(c, :label) && !isnothing(c.label[]), p.plots)
end
