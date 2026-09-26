# The names before 0.3 still work, with a deprecation warning
@testset "deprecated names" begin
    l = @test_deprecated treelayout(tree)
    @test l.names == tree_layout(tree).names
    @test @test_deprecated(hassos(tree, sos, "n1")) == has_sos(tree, sos, "n1")
    @test @test_deprecated(cladecolors(:RdYlBu)) == clade_colors(:RdYlBu)
    @test @test_deprecated(focuscolors(tree, l, "n1", :RdYlBu, :gray)) ==
        focus_colors(tree, l, "n1", :RdYlBu, :gray)
    @test @test_deprecated(speciesname("Carduelis flammea")) ==
        species_name("Carduelis flammea")
    @test @test_deprecated(selectclades(l, 2)) == select_clades(l, 2)
    @test @test_deprecated(imagegeometry(l, :fan)) == image_geometry(l, :fan)
    has = Dict("a" => true, "c" => true)
    @test @test_deprecated(imageclades(l, 2, asm, has)) == image_clades(l, 2, asm, has)

    fig, np = @test_deprecated nodepanel(asm, tree, "root", res)
    @test np isa NodePanel
    # separate figures: a node panel with an SOS map beside it overflows the stack in
    # Makie's layout, a known bug that the old names have nothing to do with
    @test @test_deprecated(nodepanel!(Figure()[1, 1], asm, tree, "root", res)) isa NodePanel
    ax, m = @test_deprecated sosmap!(Figure()[1, 1], asm, "root", res)
    @test m isa SiteMap
    et = @test_deprecated explorertree!(
        Figure()[1, 1], tree, Observable("n1"), Dict("n1" => 1.0)
    )
    @test et isa ExplorerTree

    fig, ex = @test_deprecated nodeexplorer(asm, tree, res; ordination=false)
    @test ex isa NodeExplorer
    @test @test_deprecated(nodeat(ex.treeplot, nothing, 0)) === nothing
    observer = on_node_click(n -> nothing, ex.axis, ex.treeplot)
    @test @test_deprecated(onnodeclick(n -> nothing, ex.axis, ex.treeplot)) isa
        typeof(observer)

    dir = mktempdir()
    NodivMakie.FileIO.save(joinpath(dir, "a.png"), fill(RGBAf(0.8, 0.2, 0.2, 1), 8, 8))
    fig, ax, tp = treeplot(tree; treetype=:fan, showtips=false)
    ti = @test_deprecated treeimages!(ax, tp, dir, asm)
    @test ti isa TreeImages
    @test @test_deprecated(missingimages(ti)) == missing_images(ti)
    fig, np = node_panel(asm, tree, "root", res)
    @test length(@test_deprecated(cladeimages!(np, tree, dir, asm))) == 2
end
