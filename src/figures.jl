# Ready-made figures: a site map with its colour bar, the divergence scores on the tree,
# and node panels written to a PDF.

"""
    map_figure(args...; title = "", label = "", axis = (;), kwargs...)

A [`sitemap`](@ref) of `args` in a new figure, with a colour bar labelled `label` beside
it, and `title` over the map. Other keyword arguments go to `sitemap`. Returns a
`FigureAxisPlot`, like `sitemap`.
"""
function map_figure(args...; title = "", label = "", axis = (;), kwargs...)
    fig, ax, p = sitemap(args...; axis = merge((; title), axis), kwargs...)
    Colorbar(fig[1, 2], p; label)
    return Makie.FigureAxisPlot(fig, ax, p)
end

"""
    metric_tree(tree, res; metric, nodes, kwargs...)
    metric_tree(tree, values::AbstractDict; nodes = keys(values), label = "value", kwargs...)

The tree with a marker at each of `nodes`, coloured by its value, and a colour bar: the
Makie version of Nodiv's `plot_gnd`. No other node gets a marker.

From an analysis result `res`, the value is the score `metric` (default `:rms` for a
`NodeMetrics`, `:gnd` for a `NodeAnalysis`; see Nodiv's `node_scores`), and `nodes` are by
default the divergent nodes by that score (Nodiv's `divergent_nodes`). The GND, a
proportion, gets the colour range 0 to 1. Or pass the values themselves, a Dict of node
name => value.

Keyword arguments:
- `title = ""`: the axis title
- `label`: the colour bar's label; the name of `metric`, or `"value"`
- `treetype = :fan`, `showtips = false`, `markersize = 8`, `strokewidth = 0.5`,
  `colormap = :YlOrRd`, `colorrange`, and any other keyword arguments: passed to
  [`treeplot`](@ref)
- `axis = (;)`, `figure = (;)`: attributes for the axis and the figure

Returns a `FigureAxisPlot`.
"""
function metric_tree(tree, res::Nodiv.AbstractNodeResult; metric = default_score(res),
                     nodes = automatic, label = string(metric),
                     colorrange = metric === :gnd ? (0, 1) : automatic, kwargs...)
    nodes = nodes === automatic ? defaultnodes(res, metric) : nodes
    return metric_tree(tree, node_scores(res, metric); nodes, label, colorrange, kwargs...)
end

function metric_tree(tree, values::AbstractDict; nodes = keys(values), label = "value",
                     title = "", treetype = :fan, showtips = false, markersize = 8,
                     strokewidth = 0.5, colormap = :YlOrRd, colorrange = automatic,
                     axis = (;), figure = (;), kwargs...)
    nodes = collect(nodes)
    fig, ax, p = treeplot(tree; treetype, showtips,
                          nodecolor = Dict(n => values[n] for n in nodes),
                          shownodes = nodes, markersize, strokewidth, colormap, colorrange,
                          axis = merge((; title), axis),
                          figure = merge((; size = (800, 700)), figure), kwargs...)
    Colorbar(fig[1, 2], p; label)
    return Makie.FigureAxisPlot(fig, ax, p)
end

"""
    node_panel_pdf(assemblage, tree, nodes, res, outfile; backend, kwargs...) -> outfile

Write the [`nodepanel`](@ref) of each of `nodes` to `outfile`, a PDF with one node per
page. `res` is the cached analysis result, as for `nodepanel`; nothing is recomputed. Other
keyword arguments go to `nodepanel`.

`backend` is a Makie backend that writes PDF, i.e. `CairoMakie` (NodivMakie does not load
one itself). The panel is built once and switched from node to node, each page is saved
on its own, and the pages are joined with `pdfunite`, from poppler, which must be on the
PATH (e.g. `brew install poppler`).
"""
function node_panel_pdf(assemblage, tree, nodes, res, outfile; backend, kwargs...)
    if Sys.which("pdfunite") === nothing
        throw(ErrorException("node_panel_pdf needs `pdfunite` (poppler) on the PATH, " *
                             "e.g. `brew install poppler`"))
    end
    tmp = mktempdir()
    pages = String[]
    fig, panel = nodepanel(assemblage, tree, first(nodes), res; kwargs...)
    for (i, node) in enumerate(nodes)
        panel.node[] = node
        page = joinpath(tmp, string(lpad(i, 3, '0'), ".pdf"))
        save(page, fig; backend)
        push!(pages, page)
    end
    run(`pdfunite $pages $outfile`)
    rm(tmp; recursive = true)
    @info "wrote node-panel PDF" outfile npages = length(pages)
    return outfile
end
