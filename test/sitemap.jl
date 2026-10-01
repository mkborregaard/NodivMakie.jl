@testset "sitemap" begin
    @test nsites(asm) == 12
    fig, ax, p = sitemap(asm)
    @test p isa SiteMap
    bg, hm = p.plots
    @test bg isa Heatmap && hm isa Heatmap
    @test size(p.image[]) == (4, 3)
    # richness at the grid cell (x = 1, y = 1) is species a and d
    @test p.image[][1, 1] == 2
    @test size(Makie.colorbuffer(fig)) != (0, 0)
    Colorbar(fig[1, 2], p)
    # empty sites are NaN, so not drawn
    empty = view(asm, species=["e"])
    @test count(isnan, sitemap(empty).plot.image[]) == 9
    # every site is drawn in `empty_color` under the values: transparent by default, so
    # sites with no value look like no site
    @test all(==(RGBAf(0, 0, 0, 0)), p.sites[])
    fig, ax, p = sitemap(empty; empty_color=:lightgray)
    @test count(==(to_color(:lightgray)), p.sites[]) == nsites(empty)
    @test count(==(RGBAf(0, 0, 0, 0)), p.sites[]) == length(p.sites[]) - nsites(empty)
    p.empty_color = :red
    @test count(==(to_color(:red)), p.sites[]) == nsites(empty)
    # the colour bar is the values' (the sites have no colormap)
    cb = Colorbar(fig[1, 2], p)
    @test to_colormap(cb.colormap[]) == to_colormap(p.colormap[])
    @test size(Makie.colorbuffer(fig)) != (0, 0)
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
    fig, ax, p = sitemap(pts; empty_color=:lightgray)
    bg, sc = p.plots
    @test bg isa Scatter && sc isa Scatter && bg.color[] == to_color(:lightgray)
    @test size(Makie.colorbuffer(fig)) != (0, 0)
end
