@testset "clusters" begin
    # two clusters of two and a node on its own
    D = [0 0.1 0.9 0.9 1; 0.1 0 0.9 0.9 1; 0.9 0.9 0 0.2 1; 0.9 0.9 0.2 0 1; 1 1 1 1 0]
    nodes = ["root", "n1", "n2", "n3", "a"]
    c = sos_clusters(D, nodes; simcut = 0.7)
    @test length(c.labels) == 2
    @test length(cluster_colors(25)) == 25 && cluster_colors(3) == cluster_colors(25)[1:3]

    fig = sos_cluster_heatmap(c; title = "clusters")
    @test fig isa Figure
    hm = only(filter(x -> x isa Axis && x.xticks[] isa Tuple, fig.content))
    @test hm.yticks[][2] == nodes[c.hclust.order]
    heat = only(filter(x -> x isa Heatmap, hm.scene.plots))
    @test heat[3][] ≈ (1 .- D)[c.hclust.order, c.hclust.order]
    # each cluster of more than one node outlined, twice (a black and a coloured line)
    @test count(x -> x isa Poly, hm.scene.plots) == 4
    # the dendrogram: leaves at x = 0 beside heatmap rows 1:n, merges at minus their heights
    dend = only(filter(x -> x isa Axis && x.xlabel[] == "1 - |r|", fig.content))
    linepoints(p) = p isa Union{Lines, LineSegments} ? Point2d.(p[1][]) :
                    reduce(vcat, map(linepoints, p.plots); init = Point2d[])
    pts = filter(q -> all(isfinite, q), reduce(vcat, map(linepoints, dend.scene.plots)))
    @test Set(round.(first.(pts); digits = 9)) == Set([0; -c.hclust.heights])
    @test Set(round(Int, q[2]) for q in pts if abs(q[1]) < 1e-6) == Set(1:length(nodes))
    @test size(Makie.colorbuffer(fig)) != (0, 0)

    fig, ax, p = cluster_tree(tree, c; title = "on the tree")
    @test ax.title[] == "on the tree"
    groups = filter(x -> x isa Scatter && haskey(x, :label) && !isnothing(x.label[]), p.plots)
    @test sort([g.label[] for g in groups]) == ["1", "2"]
    @test any(x -> x isa Legend, fig.content)
    @test size(Makie.colorbuffer(fig)) != (0, 0)
    # no clusters of more than one node: no markers, no legend
    fig, ax, p = cluster_tree(tree, sos_clusters(D, nodes; simcut = 0.95))
    @test !any(x -> x isa Legend, fig.content)
end
