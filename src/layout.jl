# Node coordinates for drawing a tree. Positions are those of Phylo's own plot recipe
# (`Phylo._findxy`), so figures match the Plots versions node for node, and a node
# index here is the same as a position in Phylo's plotting order.

"""
    TreeLayout

Plotting coordinates of every node in a tree, in Phylo's plotting order.

- `names`: node names
- `height`: distance from the root (the x axis of a dendrogram, the radius of a fan)
- `depth`: vertical position; tips are at `1:ntips`, internal nodes midway between
  their outermost children
- `parent`: index of each node's parent, `0` for the root
- `isleaf`: whether each node is a tip
- `index`: node name => index
"""
struct TreeLayout
    names::Vector{String}
    height::Vector{Float64}
    depth::Vector{Float64}
    parent::Vector{Int}
    isleaf::BitVector
    index::Dict{String, Int}
end

"""
    treelayout(tree)

Compute the [`TreeLayout`](@ref) of a Phylo tree.
"""
function treelayout(tree::Phylo.AbstractTree)
    h, d, n = Phylo._findxy(tree)
    names = String[getnodename(tree, x) for x in n]
    index = Dict(name => i for (i, name) in enumerate(names))
    parent = [hasinbound(tree, x) ? index[getnodename(tree, getparent(tree, x))] : 0
              for x in n]
    TreeLayout(names, Float64[h[x] for x in n], Float64[d[x] for x in n], parent,
               BitVector([isleaf(tree, x) for x in n]), index)
end

Base.length(l::TreeLayout) = length(l.names)
ntips(l::TreeLayout) = count(l.isleaf)

# Angle of a depth value in a fan tree; tips are spread over the full circle with a
# one-tip gap between the first and the last, as in Phylo.
fanangle(depth, ntip) = 2pi * depth / (ntip + 1)

polar(r, θ) = Point2d(r * cos(θ), r * sin(θ))

"""
    nodepositions(layout, treetype)

Coordinates of every node for `treetype` `:dendrogram` or `:fan`.
"""
function nodepositions(l::TreeLayout, treetype::Symbol)
    if treetype === :dendrogram
        return Point2d.(l.height, l.depth)
    elseif treetype === :fan
        nt = ntips(l)
        return [polar(r, fanangle(d, nt)) for (r, d) in zip(l.height, l.depth)]
    end
    throw(ArgumentError("Unsupported `treetype` $(repr(treetype)); valid values are `:dendrogram` or `:fan`"))
end

"""
    branchpaths(layout, treetype)

One polyline per branch, NaN-separated, for drawing with `lines`. Returns the points
and, for each point, the index of the node the branch leads to, so per-node values can
be expanded to per-vertex colours.

A dendrogram branch runs vertically along its parent's height, then horizontally to
the node. A fan branch runs along an arc at its parent's radius, then radially out.
"""
function branchpaths(l::TreeLayout, treetype::Symbol)
    pts = Point2d[]
    owner = Int[]
    nan = Point2d(NaN, NaN)
    nt = ntips(l)
    for i in eachindex(l.names)
        p = l.parent[i]
        p == 0 && continue
        hp, dp, hi, di = l.height[p], l.depth[p], l.height[i], l.depth[i]
        n0 = length(pts)
        if treetype === :dendrogram
            push!(pts, Point2d(hp, dp), Point2d(hp, di), Point2d(hi, di))
        elseif treetype === :fan
            θp, θi = fanangle(dp, nt), fanangle(di, nt)
            # same arc resolution as Phylo's `_p_circ`
            for θ in range(θp, θi; length = 1 + ceil(Int, 60abs(θi - θp)))
                push!(pts, polar(hp, θ))
            end
            push!(pts, polar(hi, θi))
        else
            nodepositions(l, treetype)  # throws the informative error
        end
        push!(pts, nan)
        append!(owner, Iterators.repeated(i, length(pts) - n0))
    end
    return pts, owner
end
