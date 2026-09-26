@testset "ordination" begin
    # root and n1 have nearly the same SOS map, as have n2 and n3
    p1, p2, e = collect(1.0:12) .- 6.5, 5 .* sin.(1:12), 0.3 .* cos.(3 .* (1:12))
    sos2 = Dict("root" => p1, "n1" => p1 .+ e, "n2" => p2, "n3" => p2 .- e)
    res2 = NodeAnalysis(
        internal, Dict("root" => 0.2, "n1" => 0.85, "n2" => 0.9, "n3" => 0.95), sos2
    )
    nodes = ["root", "n1", "n2", "n3"]
    o = sos_ordination(res2, nodes)                              # from Nodiv
    dist(a, b) = sqrt(
        sum(
            abs2,
            o.coords[:, findfirst(==(a), nodes)] .- o.coords[:, findfirst(==(b), nodes)],
        ),
    )
    @test dist("root", "n1") < dist("root", "n2")                 # similar maps plot close
    @test dist("n2", "n3") < dist("n1", "n3")
    thin = merge(sos2, Dict("n3" => [p2[1:2]; fill(NaN, 10)]))   # n3 shares two cells
    # the eigenvalues, from an ordination with more axes
    o3 = sos_ordination(o.distances, nodes; maxoutdim=3)
    fig, ax, p = eigenvalueplot(o3)
    @test p isa EigenvaluePlot && ax.xlabel[] == "MDS axis"
    @test last.(p.bars[]) == o3.eigenvalues
    @test size(Makie.colorbuffer(fig)) != (0, 0)

    fig, ax, p = ordinationplot(o; nodecolor=res2.gnd, nodelabels=true, selected="n2")
    @test p isa OrdinationPlot
    @test ax.autolimitaspect[] == 1 && ax.xlabel[] == "MDS axis 1"
    pts = p.points[]
    @test pts == Point2d.(o.coords[1, :], o.coords[2, :])
    @test p.point_colors[] == [0.2, 0.85, 0.9, 0.95]
    @test p.label_texts[] == nodes
    @test p.selected_points[] == [pts[3]]                        # a ring at n2
    p.selected = "a"                                             # not in the ordination
    @test isempty(p.selected_points[])
    p.selected = nothing
    @test isempty(p.selected_points[])
    @test nodeat(p, p.plots[1], 2) == "n1"
    @test nodeat(p, p.plots[2], 1) === nothing                   # the ring
    @test nodeat(p, p.plots[1], 5) === nothing
    @test nodeat(p, nothing, 0) === nothing
    hover(plt, i) = plt.inspector_label[](plt, i, nothing)
    @test hover(p.plots[1], 4) == "n3"
    Colorbar(fig[1, 2], p)
    @test size(Makie.colorbuffer(fig)) != (0, 0)
    fig, ax, p = ordinationplot(o; nodecolor=Dict("n1" => 1.0))
    @test count(isnan, p.point_colors[]) == 3                    # the others transparent

    # in the explorer, in place of the clade map, and linked both ways with the tree
    target = Ref{Any}((nothing, 0))
    fig, ex = nodeexplorer(
        asm, tree, res2; nodes=:all, treetype=:dendrogram, pickfn=(sc, xy, r) -> target[]
    )
    op = ex.ordination
    @test op isa OrdinationPlot && ex.ordinationaxis isa Axis
    @test ex.panel.axes[1] === nothing && length(ex.panel.colorbars) == 3
    @test op.ordination[].nodes == sort(nodes)
    @test op.selected[] == ex.panel.node[] == "n3"               # highest GND, ringed
    @test op.point_colors[] == [res2.gnd[n] for n in sort(nodes)]  # as on the tree
    @test op.colorrange[] == ex.treeplot.joint_colorrange[]
    Makie.colorbuffer(fig)
    e = events(fig)
    vp = ex.ordinationaxis.scene.viewport[]
    e.mouseposition[] = Tuple(Float64.(vp.origin .+ vp.widths ./ 2))
    target[] = (op.plots[1], findfirst(==("n1"), op.ordination[].nodes))
    e.mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.press)
    e.mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.release)
    @test ex.panel.node[] == "n1"                                # a click on a point
    @test op.selected[] == "n1"
    @test startswith(ex.status[], "n1   gnd = 0.85")
    ex.panel.node[] = "n2"                                       # the ring follows the tree
    @test op.selected[] == "n2"
    @test hover(op.plots[1], findfirst(==("n2"), op.ordination[].nodes)) ==
        "n2  (3 species)\ngnd = 0.9"
    @test size(Makie.colorbuffer(fig)) != (0, 0)
    # options for the ordination
    fig, ex = nodeexplorer(
        asm,
        tree,
        NodeAnalysis(res2.nodes, res2.gnd, thin);
        nodes=:all,
        ordinationkw=(; minoverlap=2),
    )
    @test ex.ordination.ordination[].distances[4, 3] < 1
    # the clade map instead: on request, or with fewer than three marked nodes
    fig, ex = nodeexplorer(asm, tree, res2; nodes=:all, ordination=false)
    @test ex.ordination === nothing && ex.panel.axes[1] isa Axis
    fig, ex = nodeexplorer(asm, tree, res2; nodes=["n1", "n2"])
    @test ex.ordination === nothing && ex.panel.axes[1] isa Axis
end
