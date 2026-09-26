# NodivMakie

Makie plotting for [Nodiv](https://github.com/mkborregaard/Nodiv.jl): phylogenies and
maps for node-based analysis, built with Makie so that they are interactive (clickable
nodes linked to maps, hover labels) and can be combined in complex layouts.

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

The two child-clade map titles are coloured like the clades' branches in the explorer:
the first child in the high (blue) end of the SOS colour map, the second in the low (red)
end (`titlecolors = false` to turn this off).

`nodepanel!(fig[1, 2], ...)` puts the panel into any layout position. `np.node` is an
`Observable{String}`, so anything can drive it. `clademap = false` leaves out the clade
map and keeps its cell (`np.layout[1, 1]`) free; `np.axes[1]` is then `nothing`.

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
panel beside it. In the panel, the map of the node's own clade is replaced by an
ordination of the marked nodes by the similarity of their SOS maps (see
[below](#ordination-by-sos-similarity-sosordination)). By default these are Nodiv's `divergent_nodes(res)` with its default
threshold, and the metric is `:rms` for a `NodeMetrics` or `:gnd` for a `NodeAnalysis`.
The node shown first is the most divergent one. Node markers have a thin dark outline
(`strokewidth`, `strokecolor`), so high-metric markers stay visible on the red clade. Clicking
a node marker, or the branch leading to a node, shows that node in the panel, and so does
clicking a point in the ordination; the node shown has a ring there. In the
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

Hovering shows a label (a Makie `DataInspector`, on by default with `inspector = true`).
On the tree and the ordination it gives the node under the cursor, with its number of
species and metric value; a branch gives the node it leads to. On the maps it gives the cell's value. On a bird around the tree it gives the species and
the clade it stands for, and on a bird in a map corner the species. Only the bird itself
counts: over its transparent surroundings, whatever is underneath is inspected.
Outside the explorer, call `DataInspector(fig)` after `treeplot`; the tree's labels come
from its `hoverlabel` attribute, a function of the node name. Like clicking, this needs
an interactive backend.

The ordination is computed from the cached SOS with `sos_distances`' default minimum
overlap of 3 cells. Set it for the space with `ordinationkw = (; minoverlap = 8)`. With
fewer than three marked nodes, or `ordination = false`, the panel shows the clade map.

The building blocks can be used on their own:
- `onnodeclick(f, ax, plot)` calls `f(nodename)` on a click on a tree or ordination.
- `nodeat(plot, pickedplot, index)` turns a `pick` result into a node name.
- `hassos(tree, sos, node)` tells whether a node can be shown in a panel.
- `focuscolors(tree, layout, node, sos_colormap, contextcolor)` gives the per-branch
  colours used for the selected node.
- `explorertree!(gridposition, tree, node, marked)` is the tree side of the explorer on
  its own: the marked tree, the label and the colour bar, with clicks setting `node`, an
  `Observable` of the node shown. Anything that follows `node` makes up the other panels,
  e.g. `sosmap!(gridposition, assemblage, node, res)`, the SOS map of the node shown.
  `nodeexplorer` is `explorertree!` with a `nodepanel!` beside it.

```julia
fig = Figure(size = (1600, 850))
node = Observable("Node 17672")
tr = explorertree!(fig[1, 1], tree, node, Dict(n => res_g.rms[n] for n in divergent_g);
                   label = "geo rms",
                   selectable = n -> hassos(tree, res_g.sos, n) && hassos(tree, res_e.sos, n))
sosmap!(fig[1, 2], birds_g, node, res_g; title = "Geographic SOS")
sosmap!(fig[1, 3], birds_e, node, res_e; title = "Environmental SOS")
DataInspector(fig)
```

## Ordination by SOS similarity: `sosordination`

Classical MDS of nodes by the similarity of their SOS maps, on Nodiv's `sos_distances`
(1 − |r| over the cells both nodes occupy). It uses the cached SOS in `res`.

```julia
o = sosordination(res_g, divergent_nodes(res_g); minoverlap = 8)
fig, ax, p = ordinationplot(o; nodecolor = res_g.rms, nodelabels = true)
Colorbar(fig[1, 2], p)
o.eigenvalues                                  # of the axes, largest first
```

Keyword arguments other than `maxoutdim` (the number of axes, default 2) go to
`sos_distances`: `minoverlap`, `method = :pearson | :spearman`, `overlapweight`. Set
`minoverlap` for each space, as the geographic and environmental spaces have very
different numbers of cells. With more axes (`maxoutdim = 10`) the eigenvalues show
whether two axes capture the structure. When the nodes are mostly unrelated in SOS
pattern, all distances are near 1 and the points form a ring. That is the finding, not
a failure of the method.

`ordinationplot` takes per-node `nodecolor` and `nodelabels` like `treeplot`, and
`selected` rings a node. `nodeat` and `onnodeclick` work on it as on a tree.

## Species images: `treeimages!`

Species images drawn as a ring around a fan tree, or as a column right of a dendrogram.
No images come with the package, as they are often copyrighted. You supply a directory
of image files named by species, e.g. `Carduelis_hornemanni.jpg`. The name matching
ignores case and treats spaces, hyphens and underscores alike. PNG and JPEG are read.

```julia
fig, ax, tp = treeplot(tree; treetype = :fan, showtips = false)
ti = treeimages!(ax, tp, "path/to/images", birds_g)   # range sizes from the assemblage
missingimages(ti)            # species to find images for, one per empty position

# in the explorer
fig, ex = nodeexplorer(birds_g, tree, res_g; images = "path/to/images")
```

Which species are shown:

1. **How many.** The image size sets how many fit (`imagesize`, as a fraction of the
   tree's radius for a fan, default 0.15; or pass `nimages`). For a fan, that is the
   number of discs that fit around a ring just outside the tips with a little space
   between them. That many equal slots divide the tips.
2. **Which clades.** Each image stands for a monophyletic clade and sits at its centre.
   The clades are disjoint, so no species is in two images. Each is at least `minclade`
   of a slot wide (default 0.5), and their centres are at least a slot apart, so images
   never overlap. Among the choices meeting these rules, the one with the most images is
   taken, and then the one covering the most species. So a clade too small for its own
   image joins its sister in their parent's image where that costs no image, and a clade
   wider than a slot leaves gaps beside its image. The choice is exact (see
   `selectclades`).
3. **Which species.** Each clade is shown by its species with the largest range size
   (occupied cells in the assemblage, or a Dict you pass) that has an image. Where none of
   a clade's species has an image, the position is left empty.

On the workshop's bird tree (9852 species), the defaults give 26 images around the fan
and 12 beside a dendrogram, covering 90–95% of species. `minclade = 1` makes every clade
fill at least its own slot; that gives fewer images (19 and 8).

Images are shown whole, shrunk to fit inside a disc (`shape = :square` for squares), reduced
to thumbnails, and outlined. `fit = :crop` fills the disc with the image's central part
instead, which suits photos. `whitebackground = true` makes a white background transparent,
for illustrations drawn on white, such as the Birds of the World ones; the image is then
sized so the bird itself, not its white card, fills the disc. `clip = 0.02` enlarges it
further by letting the outermost 2% of the bird (a tail tip, wingtips) be cut off at the
disc's edge; the default 0 never cuts anything. Links to missing files
are skipped, so a folder of links can point at images that are still downloading. There is no outline by default (`strokewidth = 1`
adds one). `showclades = true` marks the tips each image stands for. A dendrogram's images go in a narrow axis beside the tree, linked in
y, so they stay square however the tree axis is shaped.

### In the node panel

`cladeimages!(panel, tree, images, rangesize)` puts an image of each child clade's
widest-ranging species (with an image) in the top-right corner of its richness map. The
images follow the node shown. They are fixed in screen space, so zooming the maps leaves
them in place, and their size is `pixelsize`. `nodepanel(...; images = dir)` does this
directly. In `nodeexplorer(...; images = dir)` they are always the same size as the images
around the tree.

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
