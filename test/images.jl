@testset "species images" begin
    FileIO = NodivMakie.FileIO
    # synthetic images (none are shipped): a wide PNG, a tall JPEG, and a non-image
    dir = mktempdir()
    wide = [RGBAf(i / 300, j / 200, 0.5, 1) for i in 1:300, j in 1:200]
    FileIO.save(joinpath(dir, "a.png"), wide)
    FileIO.save(
        joinpath(dir, "B.jpg"), [Makie.RGBf(0.2, 0.4, i / 120) for i in 1:120, j in 1:80]
    )
    write(joinpath(dir, "notes.txt"), "not an image")

    @test speciesname("Carduelis_hornemanni") ==
        speciesname("carduelis  Hornemanni") ==
        speciesname("Carduelis-hornemanni") ==
        "carduelis_hornemanni"
    imgs = SpeciesImages(dir; maxpixels=64)
    @test length(imgs) == 2
    @test haskey(imgs, "a") && haskey(imgs, "b") && !haskey(imgs, "c")
    @test maximum(size(imgs["a"])) <= 64            # _thumbnail
    @test imgs["a"] === imgs["a"]                    # cached
    @test_throws KeyError imgs["c"]
    @test_throws ArgumentError SpeciesImages(joinpath(dir, "nope"))
    m = NodivMakie._marker_image(imgs["a"], :circle; fit=:crop)
    @test size(m, 1) == size(m, 2)                   # cropped to a square
    @test m[1, 1].alpha == 0 && m[end ÷ 2, end ÷ 2].alpha == 1   # disc mask
    # :pad keeps the whole image: every pixel of it is in the result, inside the disc
    h, w = size(imgs["a"])
    m = NodivMakie._marker_image(imgs["a"], :circle)
    @test size(m, 1) == size(m, 2) >= hypot(h, w)
    @test count(p -> p.alpha > 0.99, m) >= h * w - 2(h + w)
    @test size(NodivMakie._marker_image(imgs["a"], :square), 1) >= max(h, w)
    # a bird on white: once the white is transparent, the bird itself fills the disc
    bird = [
        if (i - 30)^2 / 400 + (j - 20)^2 / 100 < 1
            RGBAf(0.3, 0.2, 0.1, 1)
        else
            RGBAf(1, 1, 1, 1)
        end for i in 1:60, j in 1:40
    ]                     # a 40 x 20 oval on a 60 x 40 white card
    plain = NodivMakie._marker_image(bird, :circle)
    filled = NodivMakie._marker_image(bird, :circle; whitebackground=true)
    @test size(filled, 1) < size(plain, 1) / 1.5          # the bird is drawn much larger
    # none clipped
    @test count(p -> p.alpha > 0.99, filled) >= count(p -> p.r < 0.5, bird) - 10
    @test size(filled, 1) <= ceil(Int, 42 * 1.03) + 1     # just around the 40-pixel oval
    # clip: a larger image, with only a little of the subject cut off
    clipped = NodivMakie._marker_image(bird, :circle; whitebackground=true, clip=0.05)
    @test size(clipped, 1) < size(filled, 1)
    @test count(p -> p.alpha > 0.99, clipped) >= 0.9 * count(p -> p.r < 0.5, bird)
    # an off-white speck at the card's edge (as in JPEG illustrations) does not count
    speck = copy(bird)
    speck[1, 1] = RGBAf(0.97, 0.97, 0.97, 1)
    @test size(NodivMakie._marker_image(speck, :circle; whitebackground=true)) ==
        size(filled)
    # white made transparent
    onwhite = [i < 3 ? RGBAf(1, 1, 1, 1) : RGBAf(0.2, 0.3, 0.1, 1) for i in 1:6, j in 1:6]
    k = NodivMakie._key_white(onwhite)
    @test all(p -> p.alpha == 0, k[1:2, :]) && all(p -> p.alpha == 1, k[3:end, :])
    # links to missing files are skipped
    symlink(joinpath(dir, "nothing_here.png"), joinpath(dir, "ghost.png"))
    @test !haskey(SpeciesImages(dir), "ghost")

    # tip spans: every clade is a contiguous run of tips
    l = treelayout(tree)
    lo, hi = NodivMakie._tip_spans(l)
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
        rlo, rhi = NodivMakie._tip_spans(rl)
        ch = selectclades(rl, nslots; minclade=m, circular=true)
        w = 200 / nslots
        @test !isempty(ch)
        @test all(i -> rhi[i] - rlo[i] + 1 >= m * w - 1e-9, ch)
        @test issorted(rlo[ch])
        @test all(k -> rhi[ch[k]] < rlo[ch[k + 1]], 1:(length(ch) - 1))       # disjoint
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
        slo, shi = NodivMakie._tip_spans(sl)
        w = 9 / nslots
        cand = findall(i -> shi[i] - slo[i] + 1 >= m * w - 1e-9, eachindex(slo))
        function feasible(set)
            set = sort(set; by=i -> slo[i])
            all(
                k ->
                    shi[set[k]] < slo[set[k + 1]] &&
                    (slo[set[k + 1]] + shi[set[k + 1]] - slo[set[k]] - shi[set[k]]) / 2 >=
                    w - 1e-9,
                1:(length(set) - 1),
            )
        end
        best = maximum(0:(2 ^ length(cand) - 1)) do mask
            set = cand[[isodd(mask >> (k - 1)) for k in 1:length(cand)]]
            if feasible(set)
                (length(set), sum(i -> shi[i] - slo[i] + 1, set; init=0))
            else
                (0, 0)
            end
        end
        ch = selectclades(sl, nslots; minclade=m)
        @test (length(ch), sum(i -> shi[i] - slo[i] + 1, ch; init=0)) == best
    end
    # a lone species sister to a big clade does not pull the big clade into one image
    lone = parsenewick(
        "(x:3,(((a:1,b:1)ab:1,(c:1,d:1)cd:1)abcd:1," *
        "((e:1,f:1)ef:1,(g:1,h:1)gh:1)efgh:1)big:1)root;",
    )
    ll = treelayout(lone)
    for m in (1.0, 0.5)
        chosen = ll.names[selectclades(ll, 4; minclade=m)]
        @test !("root" in chosen) && !("big" in chosen)
        @test issubset(["abcd", "efgh"], chosen)
    end
    @test selectclades(ll, 1; minclade=1) == [ll.index["root"]]

    # representatives: largest range among the species with an image
    rs = Dict(
        "a" => 5,
        "b" => 9,
        "c" => 7,
        "d" => 1,
        "e" => 3,
        "f" => 8,
        "g" => 2,
        "h" => 4,
        "x" => 1,
    )
    has = Dict(s => true for s in ["a", "c", "d", "e", "g", "h", "x"])   # no b, no f
    cis = imageclades(ll, 4, rs, has)
    byclade = Dict(c.clade => c for c in cis)
    if haskey(byclade, "abcd")
        @test byclade["abcd"].species == "b"       # largest range overall
        @test byclade["abcd"].shown == "c"         # largest with an image
    end
    # an assemblage gives range sizes as occupancy
    rsa = NodivMakie._range_sizes(asm)
    @test rsa["a"] == 6 && rsa["e"] == 3

    # geometry: images fit with their spacing
    Random.seed!(1)
    rt = rand(Ultrametric(300))
    rl = treelayout(rt)
    g = imagegeometry(rl, :fan)
    Δ = 2pi * 300 / 301 / g.nimages
    @test 2g.radius * sin(Δ / 2) >= 1.1 * g.size - 1e-9
    @test g.size ≈ 0.15 * maximum(rl.height)
    g2 = imagegeometry(rl, :fan; nimages=20)
    @test g2.nimages == 20
    @test imagegeometry(rl, :fan; imagesize=g2.size / maximum(rl.height)).nimages >= 20
    @test imagegeometry(rl, :fan; shape=:square).nimages < g.nimages
    gd = imagegeometry(rl, :dendrogram)
    @test gd.nimages * gd.size * 1.1 <= 300 + 1e-9

    # drawing: a ring around a fan
    rdir = mktempdir()
    for t in 1:2:300          # images for every other species
        FileIO.save(
            joinpath(rdir, "tip $t.png"),
            [RGBAf(t / 300, 0.3, 1 - t / 300, 1) for i in 1:40, j in 1:40],
        )
    end
    rimgs = SpeciesImages(rdir)
    rrs = Dict("tip $t" => t for t in 1:300)
    fig, ax, tp = treeplot(rt; treetype=:fan, showtips=false)
    ti = treeimages!(ax, tp, rimgs, rrs)
    @test ti isa TreeImages
    @test length(ti.clades) >= 10
    @test all(c -> c.shown === nothing || haskey(rimgs, c.shown), ti.clades)
    @test count(p -> p isa Image, ti.plots) == count(c -> c.shown !== nothing, ti.clades)
    @test missingimages(ti) == [c.species for c in ti.clades if c.shown === nothing]
    # hovering over an image names its species and clade
    c1 = first(filter(c -> c.shown !== nothing, ti.clades))
    ip = first(filter(p -> p isa Image, ti.plots))
    @test ip.inspectable[]
    @test ip.inspector_label[](ip, (1, 1), nothing) ==
        "$(replace(c1.shown, "_" => " "))\nfor $(c1.clade) ($(length(c1.tips)) species)"
    @test ip.inspector_hover[] isa Function
    # only over the visible bird, not the transparent corners of its disc
    di = DataInspector(fig)
    Makie.colorbuffer(fig)
    (x0, x1), (y0, y1) = ip[1][].data, ip[2][].data
    function moveto(x, y)                                # data coordinates -> mouse
        px = Makie.project(ax.scene, :data, :pixel, Point3d(x, y, 0))
        events(fig).mouseposition[] = Tuple(
            Float64.(px[Vec(1, 2)] .+ ax.scene.viewport[].origin)
        )
    end
    moveto((x0 + x1) / 2, (y0 + y1) / 2)
    @test NodivMakie._over_image(ip)
    @test ip.inspector_hover[](di, ip, 1) == true         # runs against Makie's inspector
    @test di.plot.text[] == ip.inspector_label[](ip, (1, 1), nothing)
    moveto(x0 + 0.03 * (x1 - x0), y0 + 0.03 * (y1 - y0)) # a transparent corner
    @test !NodivMakie._over_image(ip)
    @test ip.inspector_hover[](di, ip, 1) == false
    @test size(Makie.colorbuffer(fig)) != (0, 0)
    # a directory path works too, with options
    fig, ax, tp = treeplot(rt; treetype=:fan, showtips=false)
    ti = treeimages!(ax, tp, rdir, rrs; nimages=12, shape=:square, showclades=true)
    @test length(ti.clades) <= 12
    @test size(Makie.colorbuffer(fig)) != (0, 0)

    # a column beside a dendrogram, square on screen
    fig, ax, tp = treeplot(rt; showtips=false)
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
        FileIO.save(
            joinpath(edir, "$sp.png"), [RGBAf(0.8, 0.2, 0.2, 1) for i in 1:20, j in 1:20]
        )
    end
    fig, ex = nodeexplorer(asm, tree, res; images=edir)
    @test ex.images isa TreeImages
    @test size(Makie.colorbuffer(fig)) != (0, 0)
    @test !any(p -> p isa Lines, ex.images.plots)          # no outline by default
    # the child-clade maps get images the same size as the tree's
    pximgs = [
        only(
            filter(p -> p isa Image && p.space[] === :pixel, ex.panel.axes[k].scene.plots)
        ) for k in (3, 4)
    ]
    w = ex.images.pixelsize[]
    @test w > 10
    for p in pximgs
        @test p.visible[] == false || -(reverse(p[1][].data)...) ≈ w
    end

    # cladeimages!: each child clade's largest-range species with an image
    eimgs = SpeciesImages(edir)
    ranges = NodivMakie._range_sizes(asm)
    best(clade) = (
        sps=filter(s -> haskey(eimgs, s), nodespecies(tree, clade));
        isempty(sps) ? nothing : first(sort(sps; by=s -> (-ranges[s], s)))
    )
    fig, np = nodepanel(asm, tree, "root", res)
    plots = cladeimages!(np, tree, eimgs, asm; pixelsize=50)
    Makie.colorbuffer(fig)
    function check(node)
        np.node[] = node
        for (k, p) in enumerate(plots)
            sp = best(getnodename(tree, getchildren(tree, node)[k]))
            @test p.visible[] == (sp !== nothing)
            sp === nothing || @test p[3][] == NodivMakie._marker_image(eimgs[sp], :circle)
        end
    end
    check("root")
    check("n2")
    check("n1")
    check("n3")
    # hover text follows the node shown
    np.node[] = "root"
    k1 = getnodename(tree, getchildren(tree, "root")[1])
    @test plots[1].inspector_label[](plots[1], (1, 1), nothing) == best(k1)
    np.node[] = "n2"
    k1 = getnodename(tree, getchildren(tree, "n2")[1])
    @test plots[1].inspector_label[](plots[1], (1, 1), nothing) == best(k1)
    # in the maps' screen space too: over the bird, not its corners
    Makie.colorbuffer(fig)
    (x0, x1), (y0, y1) = plots[1][1][].data, plots[1][2][].data
    o = np.axes[3].scene.viewport[].origin
    events(fig).mouseposition[] = (o[1] + (x0 + x1) / 2, o[2] + (y0 + y1) / 2)
    @test NodivMakie._over_image(plots[1])
    events(fig).mouseposition[] = (o[1] + x0 + 1, o[2] + y0 + 1)
    @test !NodivMakie._over_image(plots[1])
    @test best("n1") == "a"                          # b has no image
    vp = np.axes[3].scene.viewport[]
    (x0, x1), (y0, y1) = plots[1][1][].data, plots[1][2][].data
    @test x1 ≈ vp.widths[1] - 6 && y1 ≈ vp.widths[2] - 6   # in the top-right corner
    @test x1 - x0 ≈ 50 && y1 - y0 ≈ 50
    # via nodepanel's keyword
    fig, np = nodepanel(asm, tree, "root", res; images=edir, imageoptions=(; pixelsize=40))
    @test count(p -> p isa Image && p.space[] === :pixel, np.axes[4].scene.plots) == 1
    @test size(Makie.colorbuffer(fig)) != (0, 0)
    # range sizes other than the assemblage's: e, not c, stands for n2 (c, d, e)
    ranked = Dict("a" => 1, "c" => 1, "e" => 10)
    pximage(ax) = only(filter(p -> p isa Image && p.space[] === :pixel, ax.scene.plots))
    fig, np = nodepanel(
        asm, tree, "root", res; images=edir, imageoptions=(; rangesize=ranked)
    )
    @test pximage(np.axes[3])[3][] == NodivMakie._marker_image(eimgs["e"], :circle)
    fig, ex = nodeexplorer(
        asm,
        tree,
        res;
        nodes=:all,
        node="root",
        images=edir,
        imageoptions=(; rangesize=ranked),
    )
    @test pximage(ex.panel.axes[3])[3][] == NodivMakie._marker_image(eimgs["e"], :circle)
    @test "e" in [c.shown for c in ex.images.clades]
    fig, ex = nodeexplorer(asm, tree, res)
    @test ex.images === nothing
end
