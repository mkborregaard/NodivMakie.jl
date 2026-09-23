# NodivMakie

Makie plotting for [Nodiv](https://github.com/mkborregaard/Nodiv.jl): phylogenies and
maps for node-based analysis, built with Makie so that they can become interactive (clickable
nodes linked to maps) and be combined in complex layouts. This is a private package for
personal use.

`using NodivMakie` loads and re-exports Makie, Phylo, SpatialEcology and Nodiv. You pick
the backend yourself: `using CairoMakie` for files, `using GLMakie` for interactive windows.
One name clashes between them: write `SpatialEcology.boundingbox` or `Makie.boundingbox`.

## Trees: `treeplot`

A Makie recipe that does what Phylo's Plots recipe does.

```julia
using NodivMakie, CairoMakie

tree = parsenewick("((a:1,b:1)n1:1,(c:1.5,(d:1,e:1)n3:0.5)n2:1)root;")

treeplot(tree)                                   # dendrogram with tip names
treeplot(tree; treetype = :fan)

# GND on the tree, like Nodiv's `plot_gnd(tree, res)`: nodes missing from the Dict
# get NaN, which is transparent, so only the analysed nodes show up
fig, ax, p = treeplot(tree; treetype = :fan, showtips = false,
                      nodecolor = res.gnd, colormap = :YlOrRd, colorrange = (0, 1))
Colorbar(fig[1, 2], p)

# categorical node groups get a legend
fig, ax, p = treeplot(tree; nodegroup = Dict("n1" => "west", "n2" => "east"))
axislegend(ax)

# into an existing layout
fig = Figure()
treeplot!(Axis(fig[1, 1]), tree; branchcolor = Dict("n1" => 1.0, "n2" => 2.0))
```

How Phylo's Plots attributes translate:

| Phylo (Plots)                | NodivMakie                                   |
|------------------------------|----------------------------------------------|
| `treetype`                   | `treetype` (`:dendrogram`, `:fan`)           |
| `showtips`, `tipfont`        | `showtips`, `tipfontsize`, `tipcolor`, `tipfont` |
| `linecolor`, `line_z`        | `branchcolor`                                |
| `marker_z`                   | `nodecolor`                                  |
| `markersize`, `markershape`… | `markersize`, `marker`, `strokewidth`, `strokecolor` |
| `marker_group`               | `nodegroup` (+ `groupcolors`)                |
| `series_annotations`         | `nodelabels` (`true` for node names)         |
| `clims`, `color`             | `colorrange`, `colormap` (shared by branches and nodes) |

Per-node values (`branchcolor`, `nodecolor`, `markersize`, `nodegroup`, `nodelabels`) can be:

- a `Dict` of node name => value;
- a vector in `treelayout(tree).names` order, of length nnodes, ninternal or ntips;
- a `String` naming a Phylo node-data field;
- a function of the node name.

A `Symbol` is always a colour, never a node-data field. That differs from Phylo, which
accepts both.

Node markers are drawn when `nodecolor`, `markersize` or `nodegroup` is set. `shownodes`
(`:internal`, `:all`, `:tips` or a vector of names) picks the nodes. The default follows
Phylo's rule: every node if a spec covers every node, otherwise the internal nodes.

When `treeplot` creates the axis, it hides the frame, ticks and grid, and fixes the aspect
ratio for fans. When you plot into an existing `Axis`, call `hidedecorations!(ax)` and
`hidespines!(ax)` yourself.

Differences from Phylo, all visual only:

- Fan tip labels on the left half are flipped so they read left to right.
- Tip-label gaps are in pixels (`tipoffset`), not in data units.

## Towards interactivity

The plot object exposes the computed geometry. `p.tree_layout[]` is a `TreeLayout`
(node names, heights, depths, parent indices). `p.node_points[]` holds the node coordinates,
and `p.shown[]` the layout indices of the nodes with markers. So a pick on the marker
scatter maps back to a node name through
`p.tree_layout[].names[p.shown[][idx]]`.

## Testing

```julia
# with TestEnv, in the test environment
include("test/runtests.jl")
```
