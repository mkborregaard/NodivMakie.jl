using NodivMakie
using CairoMakie
using Test

const NEWICK = "((a:1,b:1)n1:1,(c:1.5,(d:1,e:1)n3:0.5)n2:1)root;"

child(p, T) = filter(c -> c isa T, p.plots)
branchlines(p) = only(child(p, Lines))
markers(p) = child(p, Scatter)[2:end]   # the first scatter is the invisible padding

@testset "NodivMakie.jl" begin
    tree = parsenewick(NEWICK)

    @testset "layout matches Phylo" begin
        l = treelayout(tree)
        h, d, n = Phylo._findxy(tree)
        @test l.names == [getnodename(tree, x) for x in n]
        @test l.height == [h[x] for x in n]
        @test l.depth == [d[x] for x in n]
        @test count(==(0), l.parent) == 1
        @test l.names[findfirst(==(0), l.parent)] == "root"
        @test l.names[l.parent[l.index["d"]]] == "n3"
        @test count(l.isleaf) == 5
    end

    @testset "dendrogram and fan render" begin
        for tt in (:dendrogram, :fan)
            fig, ax, p = treeplot(tree; treetype = tt)
            @test p isa TreePlot
            @test size(Makie.colorbuffer(fig)) != (0, 0)
            # one branch per non-root node, NaN-separated
            @test count(q -> isnan(q[1]), p.branch_points[]) == 4 + 4
            @test all(iszero, only(child(p, Scatter)).markersize[])   # no node markers by default
        end
        @test_throws ArgumentError treeplot(tree; treetype = :radial)
    end

    @testset "branch colours" begin
        vals = Dict("a" => 1.0, "n1" => 2.0, "n2" => 3.0)
        fig, ax, p = treeplot(tree; branchcolor = vals)
        c = branchlines(p).color[]
        @test length(c) == length(p.branch_points[])
        @test extrema(filter(isfinite, c)) == (1.0, 3.0)
        @test any(isnan, c)                             # branches not in the Dict
        @test p.joint_colorrange[] == (1.0, 3.0)
    end

    @testset "node markers" begin
        l = treelayout(tree)
        # scalar markersize -> internal nodes only (Phylo's rule)
        fig, ax, p = treeplot(tree; markersize = 10)
        @test sort(l.names[p.shown[]]) == ["n1", "n2", "n3", "root"]
        # Dict -> all nodes, the missing ones NaN (transparent)
        gnd = Dict("n1" => 0.2, "n3" => 0.9)
        fig, ax, p = treeplot(tree; nodecolor = gnd, treetype = :fan)
        @test length(p.shown[]) == 9
        mc = only(markers(p)).color[]
        @test count(isfinite, mc) == 2
        @test mc[findfirst(==(l.index["n3"]), p.shown[])] ≈ 0.9
        cb = Colorbar(fig[1, 2], p)
        @test cb isa Colorbar
        # explicit node list and a function spec
        fig, ax, p = treeplot(tree; nodecolor = n -> length(n), shownodes = ["a", "n2"])
        @test only(markers(p)).color[] ≈ [1.0, 2.0]
    end

    @testset "node groups and labels" begin
        grp = Dict("n1" => "west", "n2" => "east", "n3" => "east")
        fig, ax, p = treeplot(tree; nodegroup = grp, nodelabels = true)
        gs = markers(p)
        @test [s.label[] for s in gs] == ["east", "west"]
        @test length(gs[1][1][]) == 2
        leg = axislegend(ax)
        @test leg isa Legend
        labels = only(filter(c -> c isa Makie.Text && !isempty(c.text[]) &&
                                  "n1" in c.text[], p.plots))
        # `true` names the shown nodes; the tips are named by the tip labels already
        @test sort(labels.text[]) == ["n1", "n2", "n3", "root"]
    end

    @testset "attributes update" begin
        fig, ax, p = treeplot(tree)
        before = copy(p.node_points[])
        p.treetype = :fan
        @test p.node_points[] != before
        @test size(Makie.colorbuffer(fig)) != (0, 0)
    end

    # Five species on a 4 x 3 grid; `sitecolumns` = one column per site
    grid = [(x, y) for y in 1:3 for x in 1:4]
    occ = zeros(Int, 5, 12)
    occ[1, 1:6] .= 1; occ[2, 4:9] .= 1; occ[3, 7:12] .= 1; occ[4, [1, 5, 9]] .= 1; occ[5, 10:12] .= 1
    sites = ["s$i" for i in 1:12]
    asm = Assemblage(occ, Float64[first.(grid) last.(grid)], sites, ["a", "b", "c", "d", "e"])
    sos = Dict(n => collect(range(-8, 8; length = 12)) .* k
               for (k, n) in enumerate(["root", "n1", "n2", "n3"]))
    res = NodeAnalysis(Dict("root" => 0.2, "n1" => 0.4, "n2" => 0.9, "n3" => 0.5), sos)

    @testset "sitemap" begin
        @test nsites(asm) == 12
        fig, ax, p = sitemap(asm)
        @test p isa SiteMap
        hm = only(p.plots)
        @test hm isa Heatmap
        @test size(p.image[]) == (4, 3)
        # richness at the grid cell (x = 1, y = 1) is species a and d
        @test p.image[][1, 1] == 2
        @test size(Makie.colorbuffer(fig)) != (0, 0)
        Colorbar(fig[1, 2], p)
        # empty sites are NaN, so not drawn
        empty = view(asm, species = ["e"])
        @test count(isnan, sitemap(empty).plot.image[]) == 9
        # a vector, a function and reactive updates
        fig, ax, p = sitemap(collect(1.0:12), asm)
        @test p.image[][4, 3] == 12
        p[1] = collect(12.0:-1:1)
        @test p.image[][4, 3] == 1
        @test sitemap(richness, asm).plot.image[][1, 1] == 2
        # point sites become a scatter
        pts = Assemblage(occ, [rand(12) rand(12)] .* 10, sites, ["a", "b", "c", "d", "e"])
        fig, ax, p = sitemap(pts)
        @test only(p.plots) isa Scatter
        @test size(Makie.colorbuffer(fig)) != (0, 0)
    end

    @testset "nodepanel" begin
        # fast clade richness == Nodiv's richness(get_clade(...)), for every node
        cr = NodivMakie.CladeRichness(asm)
        for n in treelayout(tree).names
            @test cr(tree, n) == richness(get_clade(asm, tree, n))
        end
        fig, np = nodepanel(asm, tree, "root", res)
        @test np isa NodePanel
        # children in Phylo's order, as in Nodiv's plot_node
        kids(n) = [getnodename(tree, c) for c in getchildren(tree, n)]
        @test [ax.title[] for ax in np.axes] == ["root", "SOS", kids("root")...]
        @test np.maps[2].colorrange[] == (-8, 8)
        @test length(np.colorbars) == 4
        @test size(Makie.colorbuffer(fig)) != (0, 0)
        # the parent map is the clade's richness, children only their species
        @test np.maps[1].image[][1, 1] == 2
        @test np.maps[3].image[][1, 1] == 1        # n2 = c, d, e: only d at (1, 1)
        @test np.maps[4].image[][1, 1] == 1        # n1 = a, b: only a at (1, 1)
        # a new node redraws in place
        np.node[] = "n2"
        @test [ax.title[] for ax in np.axes] == ["n2", "SOS", kids("n2")...]
        @test np.maps[2].values[] == sos["n2"]
        @test size(Makie.colorbuffer(fig)) != (0, 0)
        # a node without SOS leaves the panel as it was
        @test_logs (:warn,) (np.node[] = "a")
        @test np.axes[1].title[] == "n2"
        @test_throws ArgumentError nodepanel(asm, tree, "a", res)
    end

    @testset "nodeexplorer" begin
        fig, ex = nodeexplorer(asm, tree, res; treetype = :dendrogram)
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
        @test nodeat(tp, branchlines(tp), findfirst(==(l.index["n1"]), tp.branch_owner[])) == "n1"
        @test nodeat(tp, child(tp, Scatter)[1], 1) === nothing   # padding
        @test nodeat(tp, nothing, 0) === nothing
        ex.panel.node[] = "n1"
        @test ex.status[] == "n1   gnd = 0.4"
        # the node shown: first child's clade in the high end of the SOS colours (positive
        # SOS = first child over-represented), second child's in the low end, rest grey
        cmap = Makie.to_colormap(:RdYlBu)
        branchcolor(n) = branchlines(tp).color[][findfirst(==(l.index[n]), tp.branch_owner[])]
        c1, c2 = [getnodename(tree, c) for c in getchildren(tree, "n1")]
        @test branchcolor(c1) == last(cmap)
        @test branchcolor(c2) == first(cmap)
        @test branchcolor("n2") == branchcolor("n3") == branchcolor("c") == to_color(:gray75)
        ex.panel.node[] = "n2"                      # follows the node shown, descendants too
        k1, k2 = [getnodename(tree, c) for c in getchildren(tree, "n2")]
        @test all(n -> branchcolor(n) == last(cmap), [k1; [getnodename(tree, d) for d in getdescendants(tree, k1)]])
        @test all(n -> branchcolor(n) == first(cmap), [k2; [getnodename(tree, d) for d in getdescendants(tree, k2)]])
        @test branchcolor("n1") == branchcolor("a") == to_color(:gray75)
        @test all(==(14), markers(tp)[1].markersize[] .* 1)   # larger node markers
        @test !any(p -> p isa Scatter && !(p.parent isa TreePlot), ex.axis.scene.plots)  # no ring
        # a custom SOS colour map carries over to the tree
        fig3, ex3 = nodeexplorer(asm, tree, res; panel = (; sos_colormap = Reverse(:RdBu)))
        l3, tp3 = ex3.treeplot.tree_layout[], ex3.treeplot
        first_child = getnodename(tree, getchildren(tree, "n2")[1])
        col = branchlines(tp3).color[][findfirst(==(l3.index[first_child]), tp3.branch_owner[])]
        @test col == last(Makie.to_colormap(Reverse(:RdBu)))
        # clicks: real Makie mouse events through the explorer's handler, with a
        # stand-in for the backend's pick (CairoMakie cannot pick)
        target = Ref{Any}((nothing, 0))
        fig, ex = nodeexplorer(asm, tree, res; treetype = :dendrogram,
                               pickfn = (sc, xy, r) -> target[])
        Makie.colorbuffer(fig)                      # lay out, so the axis has a viewport
        tp, vp, e = ex.treeplot, ex.axis.scene.viewport[], events(fig)
        press(b = Mouse.left) = (e.mousebutton[] = Makie.MouseButtonEvent(b, Mouse.press);
                                 e.mousebutton[] = Makie.MouseButtonEvent(b, Mouse.release))
        inside = Tuple(Float64.(vp.origin .+ vp.widths ./ 2))
        e.mouseposition[] = inside
        target[] = (markers(tp)[1], findfirst(==(l.index["n3"]), tp.shown[]))
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
        press(Mouse.right)                          # other buttons are ignored
        e.mouseposition[] = (1.0, 1.0); press()     # so are clicks outside the axis
        @test ex.panel.node[] == "n1"
        e.mouseposition[] = inside; target[] = (nothing, 0); press()   # and empty space
        @test ex.panel.node[] == "n1"
        @test size(Makie.colorbuffer(fig)) != (0, 0)

        # only some nodes, another metric
        fig, ex = nodeexplorer(asm, tree, res; nodes = ["n1", "n3"], metric = Dict("n1" => 1.0, "n3" => 2.0))
        @test ex.panel.node[] == "n3"
        @test count(isfinite, markers(ex.treeplot)[1].color[]) == 2
    end
end
