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
    empty = view(asm, species=["e"])
    @test count(isnan, sitemap(empty).plot.image[]) == 9
    # a vector, a function and reactive updates
    fig, ax, p = sitemap(collect(1.0:12), asm)
    @test p.image[][4, 3] == 12
    p[1] = collect(12.0:-1:1)
    @test p.image[][4, 3] == 1
    @test sitemap(richness, asm).plot.image[][1, 1] == 2
    # a site statistic by name
    withstat = Assemblage(
        occ, Float64[first.(grid) last.(grid)], sites, ["a", "b", "c", "d", "e"]
    )
    addsitestats!(withstat, collect(1.0:12), :pc1)
    @test sitemap(:pc1, withstat).plot.image[][4, 3] == 12
    # point sites become a scatter
    pts = Assemblage(occ, [rand(12) rand(12)] .* 10, sites, ["a", "b", "c", "d", "e"])
    fig, ax, p = sitemap(pts)
    @test only(p.plots) isa Scatter
    @test size(Makie.colorbuffer(fig)) != (0, 0)
end
