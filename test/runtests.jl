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
end
