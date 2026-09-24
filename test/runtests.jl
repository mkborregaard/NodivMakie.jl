using NodivMakie
using CairoMakie
using Test
using Random

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
        fig, ex = nodeexplorer(asm, tree, res; treetype = :dendrogram, nodes = :all)
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
        blue, red = Makie.interpolated_getindex(cmap, 0.85), Makie.interpolated_getindex(cmap, 0.15)
        branchcolor(n) = branchlines(tp).color[][findfirst(==(l.index[n]), tp.branch_owner[])]
        c1, c2 = [getnodename(tree, c) for c in getchildren(tree, "n1")]
        @test branchcolor(c1) == blue
        @test branchcolor(c2) == red
        @test blue != last(cmap) && red != first(cmap)      # lighter than the end colours
        @test branchcolor("n2") == branchcolor("n3") == branchcolor("c") == to_color(:gray75)
        ex.panel.node[] = "n2"                      # follows the node shown, descendants too
        k1, k2 = [getnodename(tree, c) for c in getchildren(tree, "n2")]
        @test all(n -> branchcolor(n) == blue, [k1; [getnodename(tree, d) for d in getdescendants(tree, k1)]])
        @test all(n -> branchcolor(n) == red, [k2; [getnodename(tree, d) for d in getdescendants(tree, k2)]])
        @test branchcolor("n1") == branchcolor("a") == to_color(:gray75)
        @test all(==(14), markers(tp)[1].markersize[] .* 1)   # larger node markers
        @test markers(tp)[1].strokewidth[] == 1                # with an outline
        @test markers(tp)[1].strokecolor[] == to_color(:gray20)
        @test !any(p -> p isa Scatter && !(p.parent isa TreePlot), ex.axis.scene.plots)  # no ring
        # a custom SOS colour map carries over to the tree
        fig3, ex3 = nodeexplorer(asm, tree, res; panel = (; sos_colormap = Reverse(:RdBu)))
        l3, tp3 = ex3.treeplot.tree_layout[], ex3.treeplot
        first_child = getnodename(tree, getchildren(tree, "n2")[1])
        col = branchlines(tp3).color[][findfirst(==(l3.index[first_child]), tp3.branch_owner[])]
        @test col == Makie.interpolated_getindex(Makie.to_colormap(Reverse(:RdBu)), 0.85)
        # focusinset = 0 gives the end colours
        fig4, ex4 = nodeexplorer(asm, tree, res; focusinset = 0)
        tp4, l4 = ex4.treeplot, ex4.treeplot.tree_layout[]
        @test branchlines(tp4).color[][findfirst(==(l4.index[first_child]), tp4.branch_owner[])] == last(cmap)
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
        e.mouseposition[] = (1.0, 1.0); press()     # so are clicks outside the axis
        @test ex.panel.node[] == "n1"
        e.mouseposition[] = inside; target[] = (nothing, 0); press()   # and empty space
        @test ex.panel.node[] == "n1"
        @test size(Makie.colorbuffer(fig)) != (0, 0)

        # by default only the divergent nodes are marked (Nodiv's divergent_nodes: GND > 0.8)
        fig, ex = nodeexplorer(asm, tree, res)
        shownvals(ex) = Dict(ex.treeplot.tree_layout[].names[i] => c for (i, c) in
                             zip(ex.treeplot.shown[], markers(ex.treeplot)[1].color[]) if isfinite(c))
        @test keys(shownvals(ex)) == Set(divergent_nodes(res)) == Set(["n2"])
        # and only they get markers, so no outlines around hidden ones
        @test ex.treeplot.tree_layout[].names[ex.treeplot.shown[]] == ["n2"]
        fig, ex = nodeexplorer(asm, tree, res; nodes = :all)
        @test keys(shownvals(ex)) == Set(keys(sos))
        # a NodeMetrics: RMS-SOS by default, divergent by RMS > 1.5; :pval starts at the lowest
        nm = Nodiv.NodeMetrics(res.gnd, Dict("root" => 1.0, "n1" => 2.5, "n2" => 1.8, "n3" => 1.2),
                               Dict(n => 0.0 for n in keys(sos)), Dict(n => 0.0 for n in keys(sos)),
                               Dict("root" => 0.5, "n1" => 0.01, "n2" => 0.03, "n3" => 0.2), sos)
        fig, ex = nodeexplorer(asm, tree, nm)
        @test keys(shownvals(ex)) == Set(["n1", "n2"])
        @test ex.panel.node[] == "n1"
        @test startswith(ex.status[], "n1   rms = 2.5")
        fig, ex = nodeexplorer(asm, tree, nm; metric = :pval)
        @test keys(shownvals(ex)) == Set(["n1", "n2"])      # divergent by pval < 0.05
        @test ex.panel.node[] == "n1"                        # lowest p
        @test_throws ArgumentError nodeexplorer(asm, tree, res; nodes = ["a"])

        # only some nodes, another metric
        fig, ex = nodeexplorer(asm, tree, res; nodes = ["n1", "n3"], metric = Dict("n1" => 1.0, "n3" => 2.0))
        @test ex.panel.node[] == "n3"
        @test count(isfinite, markers(ex.treeplot)[1].color[]) == 2
    end

    @testset "species images" begin
        FileIO = NodivMakie.FileIO
        # synthetic images (none are shipped): a wide PNG, a tall JPEG, and a non-image
        dir = mktempdir()
        wide = [RGBAf(i / 300, j / 200, 0.5, 1) for i in 1:300, j in 1:200]
        FileIO.save(joinpath(dir, "a.png"), wide)
        FileIO.save(joinpath(dir, "B.jpg"), [Makie.RGBf(0.2, 0.4, i / 120) for i in 1:120, j in 1:80])
        write(joinpath(dir, "notes.txt"), "not an image")

        @test speciesname("Carduelis_hornemanni") == speciesname("carduelis  Hornemanni") ==
              speciesname("Carduelis-hornemanni") == "carduelis_hornemanni"
        imgs = SpeciesImages(dir; maxpixels = 64)
        @test length(imgs) == 2
        @test haskey(imgs, "a") && haskey(imgs, "b") && !haskey(imgs, "c")
        @test maximum(size(imgs["a"])) <= 64            # thumbnail
        @test imgs["a"] === imgs["a"]                    # cached
        @test_throws KeyError imgs["c"]
        @test_throws ArgumentError SpeciesImages(joinpath(dir, "nope"))
        m = NodivMakie.markerimage(imgs["a"], :circle)
        @test size(m, 1) == size(m, 2)                   # cropped to a square
        @test m[1, 1].alpha == 0 && m[end ÷ 2, end ÷ 2].alpha == 1   # disc mask

        # tip spans: every clade is a contiguous run of tips
        l = treelayout(tree)
        lo, hi = NodivMakie.tipspans(l)
        @test (lo[l.index["root"]], hi[l.index["root"]]) == (1, 5)
        for (i, n) in enumerate(l.names)
            @test hi[i] - lo[i] + 1 == length(nodespecies(tree, n))
        end

        # the selection rule, on random trees: disjoint clades at least minclade slots
        # wide, whose centres (the image positions) are at least a slot apart
        for seed in 1:5, nslots in (3, 10, 40), m in (1.0, 0.5)
            Random.seed!(seed)
            rt = rand(Ultrametric(200))
            rl = treelayout(rt)
            rlo, rhi = NodivMakie.tipspans(rl)
            ch = selectclades(rl, nslots; minclade = m, circular = true)
            w = 200 / nslots
            @test !isempty(ch)
            @test all(i -> rhi[i] - rlo[i] + 1 >= m * w - 1e-9, ch)
            @test issorted(rlo[ch])
            @test all(k -> rhi[ch[k]] < rlo[ch[k + 1]], 1:length(ch) - 1)       # disjoint
            ctr = [(rlo[i] + rhi[i]) / 2 for i in ch]
            @test all(>=(w - 1e-9), diff(ctr))                                   # no overlap
            length(ch) > 1 && @test 201 - (ctr[end] - ctr[1]) >= w - 1e-9       # ring closes
            m == 1 && @test length(ch) <= nslots
        end
        # the most images possible (then the most species): brute force on small trees
        for seed in 1:4, nslots in (2, 3, 5), m in (1.0, 0.5)
            Random.seed!(seed)
            st = rand(Ultrametric(9))
            sl = treelayout(st)
            slo, shi = NodivMakie.tipspans(sl)
            w = 9 / nslots
            cand = findall(i -> shi[i] - slo[i] + 1 >= m * w - 1e-9, eachindex(slo))
            function feasible(set)
                set = sort(set; by = i -> slo[i])
                all(k -> shi[set[k]] < slo[set[k + 1]] &&
                         (slo[set[k + 1]] + shi[set[k + 1]] - slo[set[k]] - shi[set[k]]) / 2 >= w - 1e-9,
                    1:length(set) - 1)
            end
            best = maximum(0:(2^length(cand) - 1)) do mask
                set = cand[[isodd(mask >> (k - 1)) for k in 1:length(cand)]]
                feasible(set) ? (length(set), sum(i -> shi[i] - slo[i] + 1, set; init = 0)) : (0, 0)
            end
            ch = selectclades(sl, nslots; minclade = m)
            @test (length(ch), sum(i -> shi[i] - slo[i] + 1, ch; init = 0)) == best
        end
        # a lone species sister to a big clade does not pull the big clade into one image
        lone = parsenewick("(x:3,(((a:1,b:1)ab:1,(c:1,d:1)cd:1)abcd:1,((e:1,f:1)ef:1,(g:1,h:1)gh:1)efgh:1)big:1)root;")
        ll = treelayout(lone)
        for m in (1.0, 0.5)
            chosen = ll.names[selectclades(ll, 4; minclade = m)]
            @test !("root" in chosen) && !("big" in chosen)
            @test issubset(["abcd", "efgh"], chosen)
        end
        @test selectclades(ll, 1; minclade = 1) == [ll.index["root"]]

        # representatives: largest range among the species with an image
        rs = Dict("a" => 5, "b" => 9, "c" => 7, "d" => 1, "e" => 3, "f" => 8, "g" => 2, "h" => 4, "x" => 1)
        has = Dict(s => true for s in ["a", "c", "d", "e", "g", "h", "x"])   # no b, no f
        cis = imageclades(ll, 4, rs, has)
        byclade = Dict(c.clade => c for c in cis)
        if haskey(byclade, "abcd")
            @test byclade["abcd"].species == "b"       # largest range overall
            @test byclade["abcd"].shown == "c"         # largest with an image
        end
        # an assemblage gives range sizes as occupancy
        rsa = NodivMakie.rangesizes(asm)
        @test rsa["a"] == 6 && rsa["e"] == 3

        # geometry: images fit with their spacing
        Random.seed!(1)
        rt = rand(Ultrametric(300))
        rl = treelayout(rt)
        g = imagegeometry(rl, :fan)
        Δ = 2pi * 300 / 301 / g.nimages
        @test 2g.radius * sin(Δ / 2) >= 1.1 * g.size - 1e-9
        @test g.size ≈ 0.15 * maximum(rl.height)
        g2 = imagegeometry(rl, :fan; nimages = 20)
        @test g2.nimages == 20
        @test imagegeometry(rl, :fan; imagesize = g2.size / maximum(rl.height)).nimages >= 20
        @test imagegeometry(rl, :fan; shape = :square).nimages < g.nimages
        gd = imagegeometry(rl, :dendrogram)
        @test gd.nimages * gd.size * 1.1 <= 300 + 1e-9

        # drawing: a ring around a fan
        rdir = mktempdir()
        for t in 1:2:300          # images for every other species
            FileIO.save(joinpath(rdir, "tip $t.png"),
                        [RGBAf(t / 300, 0.3, 1 - t / 300, 1) for i in 1:40, j in 1:40])
        end
        rimgs = SpeciesImages(rdir)
        rrs = Dict("tip $t" => t for t in 1:300)
        fig, ax, tp = treeplot(rt; treetype = :fan, showtips = false)
        ti = treeimages!(ax, tp, rimgs, rrs)
        @test ti isa TreeImages
        @test length(ti.clades) >= 10
        @test all(c -> c.shown === nothing || haskey(rimgs, c.shown), ti.clades)
        @test count(p -> p isa Image, ti.plots) == count(c -> c.shown !== nothing, ti.clades)
        @test missingimages(ti) == [c.species for c in ti.clades if c.shown === nothing]
        @test size(Makie.colorbuffer(fig)) != (0, 0)
        # a directory path works too, with options
        fig, ax, tp = treeplot(rt; treetype = :fan, showtips = false)
        ti = treeimages!(ax, tp, rdir, rrs; nimages = 12, shape = :square, showclades = true)
        @test length(ti.clades) <= 12
        @test size(Makie.colorbuffer(fig)) != (0, 0)

        # a column beside a dendrogram, square on screen
        fig, ax, tp = treeplot(rt; showtips = false)
        ti = treeimages!(ax, tp, rimgs, rrs)
        @test ti.axis !== ax
        Makie.colorbuffer(fig)
        imgax = ti.axis
        pxx = imgax.scene.viewport[].widths[1] / imgax.finallimits[].widths[1]   # px per x unit
        pxy = imgax.scene.viewport[].widths[2] / imgax.finallimits[].widths[2]   # px per tip
        @test pxx * 1 ≈ pxy * ti.geometry.size rtol = 0.02      # image width == height
        @test imgax.finallimits[].origin[2] ≈ ax.finallimits[].origin[2]      # linked in y

        # in the explorer (range sizes from the assemblage)
        edir = mktempdir()
        for sp in ["a", "c", "e"]
            FileIO.save(joinpath(edir, "$sp.png"), [RGBAf(0.8, 0.2, 0.2, 1) for i in 1:20, j in 1:20])
        end
        fig, ex = nodeexplorer(asm, tree, res; images = edir)
        @test ex.images isa TreeImages
        @test size(Makie.colorbuffer(fig)) != (0, 0)
        fig, ex = nodeexplorer(asm, tree, res)
        @test ex.images === nothing
    end
end
