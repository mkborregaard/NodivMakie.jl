@testset "nodeexplorer" begin
    fig, ex = nodeexplorer(
        asm, tree, res; treetype=:dendrogram, nodes=:all, ordination=false
    )
    @test ex isa NodeExplorer
    @test ex.panel.node[] == "n2"               # highest GND
    @test startswith(ex.status[], "n2   gnd = 0.9")
    @test size(Makie.colorbuffer(fig)) != (0, 0)
    tp = ex.treeplot
    l = tp.tree_layout[]
    # picking: a node marker gives its node, a branch the node it leads to
    mk = markers(tp)[1]
    i = findfirst(==(l.index["n3"]), tp.shown[])
    @test nodeat(tp, mk, i) == "n3"
    @test nodeat(tp, branchlines(tp), findfirst(==(l.index["n1"]), tp.branch_owner[])) ==
        "n1"
    @test nodeat(tp, child(tp, Scatter)[1], 1) === nothing   # padding
    @test nodeat(tp, nothing, 0) === nothing
    ex.panel.node[] = "n1"
    @test ex.status[] == "n1   gnd = 0.4"
    # the child clades' map titles in the same colours as their branches
    c1c, c2c = cladecolors(:RdYlBu)
    @test ex.panel.axes[3].titlecolor[] == to_color(c1c)
    @test ex.panel.axes[4].titlecolor[] == to_color(c2c)
    @test ex.panel.axes[1].titlecolor[] != to_color(c1c)
    # hover labels: the node, its species count and metric value; a branch gives its node
    @test ex.inspector isa DataInspector
    hover(plt, i) = plt.inspector_label[](plt, i, nothing)
    bl = branchlines(tp)
    @test hover(bl, findfirst(==(l.index["n1"]), tp.branch_owner[])) ==
        "n1  (2 species)\ngnd = 0.4"
    @test hover(bl, findfirst(==(l.index["a"]), tp.branch_owner[])) == "a"
    mk = markers(tp)[1]
    @test hover(mk, findfirst(==(l.index["n3"]), tp.shown[])) ==
        "n3  (2 species)\ngnd = 0.5"
    # the node shown: first child's clade in the high end of the SOS colours (positive
    # SOS = first child over-represented), second child's in the low end, rest grey
    cmap = Makie.to_colormap(:RdYlBu)
    blue, red = Makie.interpolated_getindex(cmap, 0.85),
    Makie.interpolated_getindex(cmap, 0.15)
    branchcolor(n) = branchlines(tp).color[][findfirst(==(l.index[n]), tp.branch_owner[])]
    c1, c2 = [getnodename(tree, c) for c in getchildren(tree, "n1")]
    @test branchcolor(c1) == blue
    @test branchcolor(c2) == red
    @test blue != last(cmap) && red != first(cmap)      # lighter than the end colours
    @test branchcolor("n2") == branchcolor("n3") == branchcolor("c") == to_color(:gray75)
    ex.panel.node[] = "n2"                      # follows the node shown, descendants too
    k1, k2 = [getnodename(tree, c) for c in getchildren(tree, "n2")]
    @test all(
        n -> branchcolor(n) == blue,
        [k1; [getnodename(tree, d) for d in getdescendants(tree, k1)]],
    )
    @test all(
        n -> branchcolor(n) == red,
        [k2; [getnodename(tree, d) for d in getdescendants(tree, k2)]],
    )
    @test branchcolor("n1") == branchcolor("a") == to_color(:gray75)
    @test all(==(14), markers(tp)[1].markersize[] .* 1)   # larger node markers
    @test markers(tp)[1].strokewidth[] == 1                # with an outline
    @test markers(tp)[1].strokecolor[] == to_color(:gray20)
    # no ring
    @test !any(p -> p isa Scatter && !(p.parent isa TreePlot), ex.axis.scene.plots)
    # a custom SOS colour map carries over to the tree
    fig3, ex3 = nodeexplorer(asm, tree, res; panel=(; sos_colormap=Reverse(:RdBu)))
    l3, tp3 = ex3.treeplot.tree_layout[], ex3.treeplot
    first_child = getnodename(tree, getchildren(tree, "n2")[1])
    col = branchlines(tp3).color[][findfirst(==(l3.index[first_child]), tp3.branch_owner[])]
    @test col == Makie.interpolated_getindex(Makie.to_colormap(Reverse(:RdBu)), 0.85)
    # focusinset = 0 gives the end colours
    fig4, ex4 = nodeexplorer(asm, tree, res; focusinset=0)
    tp4, l4 = ex4.treeplot, ex4.treeplot.tree_layout[]
    @test branchlines(tp4).color[][findfirst(
        ==(l4.index[first_child]), tp4.branch_owner[]
    )] == last(cmap)
    # clicks: real Makie mouse events through the explorer's handler, with a
    # stand-in for the backend's pick (CairoMakie cannot pick)
    target = Ref{Any}((nothing, 0))
    fig, ex = nodeexplorer(
        asm, tree, res; treetype=:dendrogram, pickfn=(sc, xy, r) -> target[]
    )
    Makie.colorbuffer(fig)                      # lay out, so the axis has a viewport
    tp, vp, e = ex.treeplot, ex.axis.scene.viewport[], events(fig)
    press(b=Mouse.left) = (
        e.mousebutton[]=Makie.MouseButtonEvent(b, Mouse.press);
        e.mousebutton[]=Makie.MouseButtonEvent(b, Mouse.release)
    )
    inside = Tuple(Float64.(vp.origin .+ vp.widths ./ 2))
    e.mouseposition[] = inside
    target[] = (branchlines(tp), findfirst(==(l.index["n3"]), tp.branch_owner[]))
    press()
    @test ex.panel.node[] == "n3"
    @test ex.status[] == "n3   gnd = 0.5"
    target[] = (branchlines(tp), findfirst(==(l.index["n1"]), tp.branch_owner[]))
    press()
    @test ex.panel.node[] == "n1"               # a branch selects the node it leads to
    target[] = (branchlines(tp), findfirst(==(l.index["a"]), tp.branch_owner[]))
    press()
    @test ex.panel.node[] == "n1"               # a tip has no SOS: not shown...
    @test startswith(ex.status[], "a: no SOS")  # ...but reported
    target[] = (markers(tp)[1], findfirst(==(l.index["n2"]), tp.shown[]))
    @test nodeat(tp, target[]...) == "n2"       # the one marked node
    press(Mouse.right)                          # other buttons are ignored
    e.mouseposition[] = (1.0, 1.0)
    press()     # so are clicks outside the axis
    @test ex.panel.node[] == "n1"
    e.mouseposition[] = inside
    target[] = (nothing, 0)
    press()   # and empty space
    @test ex.panel.node[] == "n1"
    @test size(Makie.colorbuffer(fig)) != (0, 0)

    # by default only the divergent nodes are marked (Nodiv's divergent_nodes: GND > 0.8)
    fig, ex = nodeexplorer(asm, tree, res)
    shownvals(ex) = Dict(
        ex.treeplot.tree_layout[].names[i] => c for
        (i, c) in zip(ex.treeplot.shown[], markers(ex.treeplot)[1].color[]) if isfinite(c)
    )
    @test keys(shownvals(ex)) == Set(divergent_nodes(res)) == Set(["n2"])
    # and only they get markers, so no outlines around hidden ones
    @test ex.treeplot.tree_layout[].names[ex.treeplot.shown[]] == ["n2"]
    fig, ex = nodeexplorer(asm, tree, res; nodes=:all)
    @test keys(shownvals(ex)) == Set(keys(sos))
    # a NodeMetrics: RMS-SOS by default, divergent by RMS > 1.5; :pval starts at the lowest
    zeros_ = Dict(n => 0.0 for n in keys(sos))
    nm = NodeMetrics(
        res.nodes,
        res.gnd,
        Dict("root" => 1.0, "n1" => 2.5, "n2" => 1.8, "n3" => 1.2),
        zeros_,
        zeros_,
        Dict("root" => 0.5, "n1" => 0.01, "n2" => 0.03, "n3" => 0.2),
        Dict(n => 1.0 for n in keys(sos)),
        sos,
    )
    fig, ex = nodeexplorer(asm, tree, nm)
    @test keys(shownvals(ex)) == Set(["n1", "n2"])
    @test ex.panel.node[] == "n1"
    @test startswith(ex.status[], "n1   rms = 2.5")
    fig, ex = nodeexplorer(asm, tree, nm; metric=:pval)
    @test keys(shownvals(ex)) == Set(["n1", "n2"])      # divergent by pval < 0.05
    @test ex.panel.node[] == "n1"                        # lowest p
    @test_throws ArgumentError nodeexplorer(asm, tree, res; nodes=["a"])

    # only some nodes, another metric
    fig, ex = nodeexplorer(
        asm, tree, res; nodes=["n1", "n3"], metric=Dict("n1" => 1.0, "n3" => 2.0)
    )
    @test ex.panel.node[] == "n3"
    @test count(isfinite, markers(ex.treeplot)[1].color[]) == 2
end

@testset "explorer building blocks" begin
    # an explorer with other panels: the tree and two SOS maps following one node
    target = Ref{Any}((nothing, 0))
    fig = Figure()
    node = Observable("root")
    other = Dict("root" => sos["n2"], "n1" => sos["n1"])     # a second "space"
    et = explorertree!(
        fig[1, 1],
        tree,
        node,
        Dict("n1" => 0.4, "n2" => 0.9);
        selectable=n -> haskey(other, n),
        unselectable="not in both",
        label="geo",
        treetype=:dendrogram,
        pickfn=(sc, xy, r) -> target[],
    )
    ax1, m1 = sosmap!(fig[1, 2], asm, node, res; title=n -> "SOS $n")
    ax2, m2 = sosmap!(fig[1, 3], asm, node, other; colorbar=false)
    @test et isa ExplorerTree && et.node === node
    @test et.status[] == "root"
    @test ax1.title[] == "SOS root" && ax2.title[] == "SOS"
    @test m1.values[] == sos["root"] && m2.values[] == sos["n2"]
    Makie.colorbuffer(fig)
    tp, l, e = et.treeplot, et.treeplot.tree_layout[], events(fig)
    vp = et.axis.scene.viewport[]
    e.mouseposition[] = Tuple(Float64.(vp.origin .+ vp.widths ./ 2))
    click(n) = (
        target[]=(branchlines(tp), findfirst(==(l.index[n]), tp.branch_owner[]));
        e.mousebutton[]=Makie.MouseButtonEvent(Mouse.left, Mouse.press);
        e.mousebutton[]=Makie.MouseButtonEvent(Mouse.left, Mouse.release)
    )
    click("n1")
    @test node[] == "n1" && et.status[] == "n1   geo = 0.4"
    @test ax1.title[] == "SOS n1" && m1.values[] == m2.values[] == sos["n1"]
    click("n2")                                  # not selectable: reported, not shown
    @test node[] == "n1" && et.status[] == "n2: not in both"
    node[] = "n2"                                # no SOS in `other`: that map stays
    @test m1.values[] == sos["n2"] && m2.values[] == sos["n1"]
    @test tp.hoverlabel[]("n2") == "n2  (3 species)\ngeo = 0.9"
    @test size(Makie.colorbuffer(fig)) != (0, 0)
    @test_throws ArgumentError sosmap!(fig[2, 1], asm, "n2", other)
    @test_throws ArgumentError explorertree!(
        fig[2, 2], tree, node, Dict("n1" => 1.0); images=SpeciesImages(mktempdir())
    )
end

@testset "linked explorers" begin
    # the same tree in a second space, where only some nodes have an SOS
    other = NodeAnalysis(internal, res.gnd, Dict("root" => sos["n2"], "n1" => sos["n1"]))
    fig1, ex1 = nodeexplorer(asm, tree, res; nodes=:all, ordination=false)
    fig2, ex2 = nodeexplorer(asm, tree, other; nodes=:all, ordination=false)
    links = link_explorers!(tree, ex1, ex2)
    @test length(links) == 2
    ex1.panel.node[] = "n1"
    @test ex2.panel.node[] == "n1"
    ex2.panel.node[] = "root"
    @test ex1.panel.node[] == "root"
    ex1.panel.node[] = "n3"                   # no SOS in the other space: it stays
    @test ex2.panel.node[] == "root"
    off.(links)
    ex1.panel.node[] = "n1"
    @test ex2.panel.node[] == "root"
end
