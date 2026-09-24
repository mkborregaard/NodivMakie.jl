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

## Maps: `sitemap`

The Makie version of EcoBase's Plots recipes. It draws one value per site as a
heatmap for gridded sites, or as a scatter for point sites.

```julia
sitemap(birds_g)                          # richness, empty cells not drawn
sitemap(res_g.sos[node], birds_g; colormap = :RdYlBu, colorrange = (-8, 8))
sitemap(:PC1, birds_g)                    # a site statistic
sitemap(occupancy, birds_g)               # f(assemblage)
```

Missing and NaN values are transparent (`nan_color`). The values are the plot's first
argument, so `p[1] = newvalues` redraws the map in place.

## The node panel: `nodepanel`

The Makie version of Nodiv's `plot_node`. It is a 2×2 grid:
- top left: the richness of the node's clade;
- top right: the node's SOS (`RdYlBu`, −8 to 8);
- bottom row: the richness of the two descendant clades.

The four map axes are linked, so zooming one zooms them all. The panel reads the
cached SOS in `res` and never recomputes the analysis.

```julia
fig, np = nodepanel(birds_g, tree, "Node 15422", res_g)
np.node[] = "Node 17672"                  # redraws in place
```

`nodepanel!(fig[1, 2], ...)` puts the panel into any layout position. `np.node` is an
`Observable{String}`, so anything can drive it.

For speed, clade richness uses a sites × species index built once per panel. It gives
exactly `richness(get_clade(...))`, but a node switch takes about 10 ms instead of about
0.8 s on the 18k-cell geographic data.

## Linked tree and maps: `nodeexplorer`

```julia
using GLMakie
fig, ex = nodeexplorer(birds_g, tree, res_g)
# or another threshold, or every node with an SOS
fig, ex = nodeexplorer(birds_g, tree, res_g; nodes = divergent_nodes(res_g; threshold = 2))
fig, ex = nodeexplorer(birds_g, tree, res_g; nodes = :all)
```

The fan tree marks the divergent nodes, coloured by the divergence metric, with the node
panel beside it. By default these are Nodiv's `divergent_nodes(res)` with its default
threshold, and the metric is `:rms` for a `NodeMetrics` or `:gnd` for a `NodeAnalysis`.
The node shown first is the most divergent one. Node markers have a thin dark outline
(`strokewidth`, `strokecolor`), so high-metric markers stay visible on the red clade. Clicking
a node marker, or the branch leading to a node, shows that node in the panel. In the
tree, the selected node's first child clade is drawn in the high (blue) end of the SOS
colour map and its second child clade in the low (red) end. Both colours are taken a
little in from the ends of the map, so they are lighter (`focusinset`, 0 for the end
colours). The rest of the tree is greyed out (`contextcolor`). This follows Nodiv's SOS, which is the first child's
richness against the null: blue cells on the SOS map are where the blue clade is
over-represented, and red cells where the red clade is. The label above the tree gives
the node's name and metric value.
Nodes without SOS, such as tips, are reported in the label and not shown. With
CairoMakie you get the static figure for the first node, which is the one with the
highest metric.

The building blocks can be used on their own:
- `onnodeclick(f, ax, treeplot)` calls `f(nodename)` on a click.
- `nodeat(treeplot, plot, index)` turns a `pick` result into a node name.
- `hassos(tree, sos, node)` tells whether a node can be shown in a panel.
- `focuscolors(tree, layout, node, sos_colormap, contextcolor)` gives the per-branch
  colours used for the selected node.

## Plot geometry

A `treeplot` keeps its computed geometry on the plot object:
- `p.tree_layout[]` is a `TreeLayout`, with node names, heights, depths and parent indices;
- `p.node_points[]` holds the node coordinates;
- `p.shown[]` holds the layout indices of the nodes with markers;
- `p.branch_owner[]` gives, for each branch vertex, the node that branch leads to.

`nodeat` is built on these.

## Testing

```julia
# with TestEnv, in the test environment
include("test/runtests.jl")
```
