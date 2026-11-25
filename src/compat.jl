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

# Hover labels can be rich text (from `rich`, e.g. with coloured parts), kept as it is;
# anything else becomes a string. `RichText`, the type `rich` returns, is not public.
_label_text(x::Makie.RichText) = x
_label_text(x) = string(x)

# The `inspector_label`, `inspector_hover` and `inspector_clear` of a plot whose hover text,
# `label(plot, index, position)`, may be rich text. The `DataInspector`'s tooltip takes only
# strings (Makie 0.24 fixes its text type at its first, empty, label), so a string label
# goes the usual way (`show_data`), and rich text goes to a second tooltip of the
# inspector's, at the cursor, hidden again when the cursor moves on. Not public:
# `show_data`, `clear_temporary_plots!` (which makes the plot the inspector's selection, so
# that the inspector does not clear it again just after) and the inspector's fields, and
# the placement rules are copied from `update_tooltip_alignment!`.
function _hover_attributes(label)
    plain(plt, i, pos) = (t=label(plt, i, pos); t isa Makie.RichText ? String(t) : t)
    function hover(inspector, plt, i)
        t = label(plt, i, nothing)
        if !(t isa Makie.RichText)
            _hide_rich_tooltip!(inspector)
            return Makie.show_data(inspector, plt, i)
        end
        inspector.selection === plt || Makie.clear_temporary_plots!(inspector, plt)
        _show_rich_tooltip!(inspector, t)
        return true
    end
    clear(inspector, _) = _hide_rich_tooltip!(inspector)
    return (; inspector_label=plain, inspector_hover=hover, inspector_clear=clear)
end

# Each `DataInspector`'s rich-text tooltip, made on first use
const _RICH_TOOLTIPS = WeakKeyDict{DataInspector,Any}()

function _hide_rich_tooltip!(inspector)
    tt = get(_RICH_TOOLTIPS, inspector, nothing)
    tt === nothing || (tt.visible[] = false)
    return nothing
end

# An empty rich text has no glyphs, which Makie 0.24 cannot lay out: it starts with " ",
# and an empty label shows no tooltip.
function _show_rich_tooltip!(inspector, text)
    isempty(String(text)) && return _hide_rich_tooltip!(inspector)
    root = inspector.root
    tt = get!(_RICH_TOOLTIPS, inspector) do
        p = Makie.tooltip!(
            root, Observable(Point2f(0)); text=Makie.rich(" "), visible=false
        )
        translate!(p, 0, 0, inspector.attributes.depth[])
        return p
    end
    inspector.plot.visible[] = false
    wx, wy = Makie.widths(Makie.viewport(root)[])
    px, py = pos = Point2f(mouseposition_px(root))
    placement = py < 0.75wy ? (:above) : (:below)
    px < 0.25wx && (placement = :right)
    px > 0.75wx && (placement = :left)
    Makie.update!(tt; arg1=pos, placement, text, visible=true)
    return nothing
end
