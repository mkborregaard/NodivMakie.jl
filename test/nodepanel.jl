@testset "node_panel" begin
    fig, np = node_panel(asm, tree, "root", res)
    @test np isa NodePanel
    @test np.sos === sos
    # children in Phylo's order, as in Nodiv's plot_node
    kids(n) = [getnodename(tree, c) for c in getchildren(tree, n)]
    @test [ax.title[] for ax in np.axes] == ["root", "SOS", kids("root")...]
    @test (np.axes[3].titlecolor[], np.axes[4].titlecolor[]) ==
        to_color.(clade_colors(:RdYlBu))
    # the clade colours: `inset` in from the ends of the colour map, interpolated linearly
    cm = Makie.to_colormap(:RdYlBu)
    at(x) = (
        t=x * (length(cm) - 1) + 1;
        i=floor(Int, t);
        cm[i] * (1 - (t - i)) + cm[i + 1] * (t - i)
    )
    @test all(clade_colors(:RdYlBu) .≈ (at(0.85), at(0.15)))
    @test all(clade_colors(:RdYlBu; inset=0.3) .≈ (at(0.7), at(0.3)))
    @test clade_colors(:RdYlBu; inset=0) == (last(cm), first(cm))
    fig2, np2 = node_panel(asm, tree, "root", res; titlecolors=false)
    @test np2.axes[3].titlecolor[] == np2.axes[1].titlecolor[]
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
    @test_throws ArgumentError node_panel(asm, tree, "a", res)
    # without the clade map: its cell is left free, the other maps keep their positions
    fig, np = node_panel(asm, tree, "root", res; clademap=false)
    @test np.axes[1] === nothing && np.maps[1] === nothing
    @test [ax.title[] for ax in np.axes[2:4]] == ["SOS", kids("root")...]
    @test length(np.colorbars) == 3
    @test isempty(contents(np.layout[1, 1]))
    np.node[] = "n2"
    @test np.axes[2].title[] == "SOS" && np.maps[2].values[] == sos["n2"]
    @test size(Makie.colorbuffer(fig)) != (0, 0)
    # an SOS map beside the panel narrows its maps to sizes where Makie's equal-scale
    # limits and automatic tick label space have no fixed point (a StackOverflowError)
    for colorbar in (true, false)
        fig = Figure()
        np = node_panel!(fig[1, 1], asm, tree, "root", res)
        ax, m = sos_map!(fig[1, 2], asm, "root", res; colorbar)
        @test size(Makie.colorbuffer(fig)) != (0, 0)
    end
    # the sites with no SOS are light grey on the SOS maps, not on the richness maps
    fig = Figure()
    np = node_panel!(fig[1, 1], asm, tree, "root", res)
    ax, m = sos_map!(fig[1, 2], asm, "root", res)
    @test np.maps[2].empty_color[] == :lightgray && m.empty_color[] == :lightgray
    @test np.maps[3].empty_color[] == :transparent
    fig, np = node_panel(asm, tree, "root", res; sos_empty_color=:red)
    @test np.maps[2].empty_color[] == :red
end
