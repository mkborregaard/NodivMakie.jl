@testset "ready-made figures" begin
    fap = map_figure(asm; title="richness", label="species")
    @test fap isa Makie.FigureAxisPlot && fap.plot isa SiteMap
    @test fap.axis.title[] == "richness"
    @test only(filter(x -> x isa Colorbar, fap.figure.content)).label[] == "species"
    @test size(Makie.colorbuffer(fap.figure)) != (0, 0)

    # scores on the tree: markers only at the nodes asked for
    fig, ax, p = metric_tree(tree, res; nodes=["n1", "n2"], title="GND")
    @test p isa TreePlot && ax.title[] == "GND"
    l = p.tree_layout[]
    @test sort(l.names[p.shown[]]) == ["n1", "n2"]
    @test p.joint_colorrange[] == (0, 1)             # GND is a proportion
    @test only(filter(x -> x isa Colorbar, fig.content)).label[] == "gnd"
    @test size(Makie.colorbuffer(fig)) != (0, 0)
    fig, ax, p = metric_tree(tree, res)              # by default the divergent nodes
    @test p.tree_layout[].names[p.shown[]] == divergent_nodes(res)
    fig, ax, p = metric_tree(tree, Dict("n1" => 2.0, "n3" => 4.0); label="rms")
    @test p.joint_colorrange[] == (2.0, 4.0)
    @test sort(p.tree_layout[].names[p.shown[]]) == ["n1", "n3"]

    # the node panels of several nodes in one PDF
    if Sys.which("pdfunite") !== nothing
        out = joinpath(mktempdir(), "panels.pdf")
        @test node_panel_pdf(asm, tree, ["root", "n2"], res, out; backend=CairoMakie) == out
        @test isfile(out) && filesize(out) > 0
    end
end
