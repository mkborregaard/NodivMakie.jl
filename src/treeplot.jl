# Makie port of Phylo's Plots recipe (Phylo/src/plot.jl). The mapping from the Plots
# attributes:
#
#   treetype              -> treetype            (:dendrogram or :fan)
#   showtips / tipfont    -> showtips / tipfontsize, tipcolor, tipfont
#   linecolor / line_z    -> branchcolor         (a colour, or per-node values)
#   marker_z              -> nodecolor           (a colour, or per-node values)
#   markersize, marker... -> markersize, marker, strokewidth, strokecolor
#   marker_group          -> nodegroup           (one labelled scatter per group)
#   series_annotations    -> nodelabels
#   clims / color         -> colorrange / colormap, shared by branches and nodes
#
# Per-node specifications (branchcolor, nodecolor, markersize, nodegroup, nodelabels)
# accept a Dict of node name => value (nodes not in it are left out / drawn with
# `nan_color`), a vector in `treelayout(tree).names` order (length nnodes, or ninternal /
# ntips for just those nodes), a String naming a Phylo node-data field, or a function of
# the node name. A Symbol is always a colour name, never a node-data field.

"""
    treeplot(tree; kwargs...)

Plot a Phylo tree as a `:dendrogram` or `:fan`, optionally with coloured branches,
node markers, node groups and node labels. See the attribute list below; per-node
values are given as a `Dict` of node name => value, a vector in
`treelayout(tree).names` order, a `String` node-data field, or a function of the
node name.

Node markers are drawn when any of `nodecolor`, `markersize` or `nodegroup` is set.
Which nodes get them is controlled by `shownodes`.
"""
@recipe TreePlot (tree,) begin
    "Tree shape: `:dendrogram` or `:fan`."
    treetype = :dendrogram

    "Branch colour: a colour, or per-node values (the branch leading to each node) mapped through `colormap`."
    branchcolor = :black
    "Colour of branches whose value is missing or NaN."
    branch_nan_color = :gray80
    "Branch line width."
    linewidth = 1
    "Branch line style."
    linestyle = nothing

    "Node marker colour: a colour, or per-node values mapped through `colormap`. Missing/NaN values use `nan_color` (transparent by default, i.e. hidden)."
    nodecolor = nothing
    "Node marker size: a number or per-node values."
    markersize = nothing
    "Node marker shape."
    marker = :circle
    strokewidth = 0
    strokecolor = :black
    """
    Nodes that get markers and labels: `:internal`, `:all`, `:tips`, or a vector of node names.
    `automatic` picks `:all` when a per-node specification covers every node (a Dict,
    function, node-data field or length-nnodes vector), `:tips` for a length-ntips vector,
    and `:internal` otherwise - the same rule as Phylo.
    """
    shownodes = automatic
    "Categorical per-node values; draws one labelled scatter per group (for `Legend`/`axislegend`). Overrides `nodecolor`."
    nodegroup = nothing
    "Colours cycled over the sorted groups of `nodegroup`."
    groupcolors = Makie.wong_colors()

    "Text at nodes: `true` for node names (tips excluded while `showtips`), or per-node strings."
    nodelabels = nothing
    nodelabelsize = 9
    nodelabelcolor = :black
    nodelabelalign = (:right, :bottom)
    "Node label offset in pixels."
    nodelabeloffset = (-3, 1)

    "Draw tip names."
    showtips = true
    tipfontsize = 7
    tipcolor = :black
    tipfont = :regular
    "Gap between a tip and its name, in pixels."
    tipoffset = 3
    """
    Extra room for tip names, as a fraction of the tree's height (dendrogram) or radius
    (fan). `automatic` uses Phylo's margins: 0.15 and 0.5 with tips, none without.
    """
    tippad = automatic

    Makie.mixin_colormap_attributes()...
    Makie.mixin_generic_plot_attributes()...
end

Makie.convert_arguments(::Type{<:TreePlot}, tree::Phylo.AbstractTree) = (tree,)

## Resolving per-node specifications -------------------------------------------------

iscolorspec(x) = x isa Union{Symbol, Colorant, Tuple{Union{Symbol, Colorant}, Real}}

# One value per node in layout order, `missing` where the spec has no value.
pernode(x::AbstractDict, tree, l) = Any[get(x, name, missing) for name in l.names]
pernode(x::AbstractString, tree, l) = Any[getnodedata(tree, name)[x] for name in l.names]
pernode(f::Function, tree, l) = Any[f(name) for name in l.names]
function pernode(x::AbstractVector, tree, l)
    n = length(x)
    n == length(l) && return Vector{Any}(x)
    ret = Vector{Any}(missing, length(l))
    if n == count(!, l.isleaf)
        ret[.!l.isleaf] .= x
    elseif n == ntips(l)
        ret[l.isleaf] .= x
    else
        throw(ArgumentError("A per-node vector must have one value per node ($(length(l))), " *
                            "internal node ($(count(!, l.isleaf))) or tip ($(ntips(l))); got $n"))
    end
    return ret
end

# Numeric values -> Float64 with NaN for missing (to go through the colormap);
# anything else is taken as one colour per node.
function colorvalues(vals, nancolor)
    if all(v -> v isa Union{Real, Missing, Nothing}, vals)
        return Float64[v isa Real ? v : NaN for v in vals]
    end
    nc = to_color(nancolor)
    return RGBAf[v isa Union{Missing, Nothing} ? nc : to_color(v) for v in vals]
end

# Does this spec give a value for every node (-> :all), tips only (-> :tips)?
function spec_coverage(x, l)
    x isa Union{AbstractDict, AbstractString, Function} && return :all
    if x isa AbstractVector
        length(x) == length(l) && return :all
        length(x) == ntips(l) && return :tips
    end
    return :internal
end

function shown_indices(l::TreeLayout, shownodes, specs...)
    if shownodes === automatic
        covs = [spec_coverage(s, l) for s in specs if !(s === nothing || s isa Bool || iscolorspec(s))]
        shownodes = :all in covs ? :all : :tips in covs ? :tips : :internal
    end
    shownodes === :all && return collect(eachindex(l.names))
    shownodes === :internal && return findall(!, l.isleaf)
    shownodes === :tips && return findall(l.isleaf)
    shownodes isa AbstractVector && return [l.index[string(n)] for n in shownodes]
    throw(ArgumentError("`shownodes` must be :internal, :all, :tips or a vector of node names"))
end

## The recipe --------------------------------------------------------------------------

function Makie.plot!(p::TreePlot)
    # fail early with a plain error, rather than from inside the compute graph
    p.treetype[] in (:dendrogram, :fan) ||
        throw(ArgumentError("Unsupported `treetype` $(repr(p.treetype[])); valid values are `:dendrogram` or `:fan`"))

    map!(treelayout, p, :tree, :tree_layout)
    map!(p, [:tree_layout, :treetype], [:branch_points, :branch_owner]) do l, tt
        return branchpaths(l, tt)
    end
    map!(nodepositions, p, [:tree_layout, :treetype], :node_points)

    map!(p, [:tree, :tree_layout, :branch_owner, :branchcolor, :branch_nan_color],
         :branch_colors) do tree, l, owner, bc, nanc
        iscolorspec(bc) && return bc
        return colorvalues(pernode(bc, tree, l), nanc)[owner]
    end

    # The nodes with markers/labels and their per-node attributes
    map!(p, [:tree, :tree_layout, :node_points, :shownodes, :nodecolor, :markersize,
             :nodegroup, :nodelabels, :nan_color, :showtips],
         [:shown, :marker_points, :marker_colors, :marker_sizes, :marker_groups,
          :label_points, :label_texts]) do tree, l, pts, shownodes, nc, ms, grp, lab, nanc, st
        shown = shown_indices(l, shownodes, nc, ms, grp, lab)
        colors = nc === nothing ? :black :
                 iscolorspec(nc) ? nc : colorvalues(pernode(nc, tree, l)[shown], nanc)
        sizes = ms === nothing ? 8.0 : ms isa Real ? Float64(ms) :
                Float64[v isa Real ? v : 0.0 for v in pernode(ms, tree, l)[shown]]
        groups = grp === nothing ? nothing : pernode(grp, tree, l)[shown]
        texts = lab === nothing || lab === false ? fill("", length(shown)) :
                # tips already carry their names when showtips is on
                lab === true ? [st && l.isleaf[i] ? "" : l.names[i] for i in shown] :
                String[v isa Union{Missing, Nothing} ? "" : string(v)
                       for v in pernode(lab, tree, l)[shown]]
        haslabel = .!isempty.(texts)
        return shown, pts[shown], colors, sizes, groups, pts[shown][haslabel], texts[haslabel]
    end

    map!(p, [:tree_layout, :treetype, :tipoffset],
         [:tip_points, :tip_texts, :tip_rotations, :tip_aligns, :tip_offsets]) do l, tt, off
        tips = findall(l.isleaf)
        if tt === :fan
            θs = fanangle.(l.depth[tips], ntips(l))
            # names on the left half are flipped so they read left to right
            flip = cos.(θs) .< 0
            rot = [f ? θ + pi : θ for (θ, f) in zip(θs, flip)]
            align = [f ? (:right, :center) : (:left, :center) for f in flip]
            offs = [Vec2f(off * cos(θ), off * sin(θ)) for θ in θs]
            pts = [polar(r, θ) for (r, θ) in zip(l.height[tips], θs)]
        else
            rot = zeros(length(tips))
            align = fill((:left, :center), length(tips))
            offs = fill(Vec2f(off, 0), length(tips))
            pts = Point2d.(l.height[tips], l.depth[tips])
        end
        return pts, l.names[tips], rot, align, offs
    end

    # Text does not count towards axis limits; two invisible points widen them to fit
    # the tip names, like the xlims/ylims set by Phylo's recipe.
    map!(p, [:tree_layout, :treetype, :showtips, :tippad], :pad_points) do l, tt, st, pad
        pad = pad === automatic ? (st ? (tt === :fan ? 0.5 : 0.15) : 0.0) : pad
        H = maximum(l.height)
        if tt === :fan
            R = H * (1 + pad)
            return [Point2d(-R, -R), Point2d(R, R)]
        end
        return [Point2d(-0.05H, minimum(l.depth)), Point2d(H * (1 + pad), maximum(l.depth))]
    end

    # Branches and nodes share one colour scale
    map!(p, [:branch_colors, :marker_colors, :colorrange], :joint_colorrange) do bc, mc, cr
        cr === automatic || return cr
        vals = Float64[]
        for c in (bc, mc)
            c isa AbstractVector{<:Real} && append!(vals, filter(isfinite, c))
        end
        isempty(vals) && return (0.0, 1.0)
        lo, hi = extrema(vals)
        return lo == hi ? (lo - 0.5, hi + 0.5) : (lo, hi)
    end

    cmap = (colormap = p.colormap, colorscale = p.colorscale, colorrange = p.joint_colorrange,
            lowclip = p.lowclip, highclip = p.highclip, alpha = p.alpha)

    scatter!(p, p.pad_points; color = :transparent, markersize = 0, inspectable = false)

    lines!(p, p.branch_points; color = p.branch_colors, linewidth = p.linewidth,
           linestyle = p.linestyle, nan_color = p.branch_nan_color, cmap...)

    drawmarkers = any(!isnothing, (p.nodecolor[], p.markersize[], p.nodegroup[]))
    if drawmarkers && p.nodegroup[] === nothing
        scatter!(p, p.marker_points; color = p.marker_colors, markersize = p.marker_sizes,
                 marker = p.marker, strokewidth = p.strokewidth,
                 strokecolor = p.strokecolor, nan_color = p.nan_color, cmap...)
    elseif drawmarkers
        # Groups are fixed when the plot is created, as they set the number of scatters
        groups = sort(unique(skipmissing(p.marker_groups[])))
        for (k, g) in enumerate(groups)
            pts, sizes = Symbol(:group_points_, k), Symbol(:group_sizes_, k)
            map!(p, [:marker_points, :marker_sizes, :marker_groups], [pts, sizes]) do mp, ms, mg
                keep = isequal.(mg, g)
                return mp[keep], ms isa AbstractVector ? ms[keep] : ms
            end
            color = p.groupcolors[][mod1(k, length(p.groupcolors[]))]
            scatter!(p, getproperty(p, pts); color, markersize = getproperty(p, sizes),
                     marker = p.marker, strokewidth = p.strokewidth,
                     strokecolor = p.strokecolor, label = string(g))
        end
    end

    text!(p, p.label_points; text = p.label_texts, fontsize = p.nodelabelsize,
          color = p.nodelabelcolor, align = p.nodelabelalign, offset = p.nodelabeloffset,
          inspectable = false)

    text!(p, p.tip_points; text = p.tip_texts, rotation = p.tip_rotations,
          align = p.tip_aligns, offset = p.tip_offsets, fontsize = p.tipfontsize,
          color = p.tipcolor, font = p.tipfont, visible = p.showtips, inspectable = false)
    return p
end

# Axis defaults, like Phylo's `framestyle = :none, grid = false` and the fan's
# `aspect_ratio = 1`. Only applied when `treeplot` creates the axis.
function Makie.preferred_axis_attributes(::Type{Axis}, p::TreePlot)
    return (; xgridvisible = false, ygridvisible = false,
            xticksvisible = false, yticksvisible = false,
            xticklabelsvisible = false, yticklabelsvisible = false,
            leftspinevisible = false, rightspinevisible = false,
            topspinevisible = false, bottomspinevisible = false,
            autolimitaspect = p.treetype[] === :fan ? 1 : nothing)
end

# `Colorbar(fig[1, 2], p)`: the node colours if they are numeric, else the branches'.
function Makie.extract_colormap(p::TreePlot)
    for child in Iterators.reverse(p.plots)
        child isa Union{Scatter, Lines} || continue
        cm = Makie.extract_colormap(child)
        cm isa Makie.ColorMapping && return cm
    end
    return nothing
end

# Legend entries are the node groups (if any), not the tree itself.
Makie.get_plots(p::TreePlot) = filter(c -> c isa Scatter && haskey(c, :label) &&
                                           !isnothing(c.label[]), p.plots)
