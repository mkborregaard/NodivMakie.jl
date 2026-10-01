# Makie port of the EcoBase/SpatialEcology Plots recipes that map a per-site variable
# over an assemblage's locations: a heatmap for gridded sites, a scatter for points.

"""
    sitemap(values, assemblage)
    sitemap(assemblage)
    sitemap(f::Function, assemblage)
    sitemap(stat::Symbol, assemblage)
    sitemap(values, locations)

Map one value per site over the sites of an assemblage. Gridded sites are drawn as a
heatmap, point sites as a scatter. `sitemap(assemblage)` maps species richness with
empty sites left out (as EcoBase's recipe does), `sitemap(f, assemblage)` maps
`f(assemblage)`, and `sitemap(stat, assemblage)` maps the site statistic column `stat`.
Missing and NaN values are drawn with `nan_color`, transparent by default. Under the
values every site is drawn in `empty_color`, so with `nan_color` transparent the sites
with no value show in `empty_color` and stand out from the cells that are not sites.
"""
@recipe SiteMap (values, locations) begin
    "Marker size for point sites."
    markersize = 6
    "Marker shape for point sites."
    marker = :rect
    "Colour of every site, drawn under the values: shows the sites with no value."
    empty_color = :transparent
    Makie.mixin_colormap_attributes()...
    Makie.mixin_generic_plot_attributes()...
end

_site_values(v) = Float64[ismissing(x) || x === nothing ? NaN : Float64(x) for x in v]

# Richness with empty sites as NaN, so they are not drawn
function _occupied(richnesses)
    r = _site_values(richnesses)
    r[r .== 0] .= NaN
    return r
end

_site_locations(asm::EcoBase.AbstractAssemblage) = getcoords(places(asm))

function Makie.convert_arguments(
    ::Type{<:SiteMap}, v::AbstractVector, asm::EcoBase.AbstractAssemblage
)
    return (_site_values(v), _site_locations(asm))
end
function Makie.convert_arguments(
    ::Type{<:SiteMap}, v::AbstractVector, locs::EcoBase.AbstractLocationData
)
    return (_site_values(v), locs)
end
function Makie.convert_arguments(T::Type{<:SiteMap}, asm::EcoBase.AbstractAssemblage)
    return (_occupied(richness(asm)), _site_locations(asm))
end
function Makie.convert_arguments(
    T::Type{<:SiteMap}, f::Function, asm::EcoBase.AbstractAssemblage
)
    return Makie.convert_arguments(T, f(asm), asm)
end
function Makie.convert_arguments(
    T::Type{<:SiteMap}, stat::Symbol, asm::EcoBase.AbstractAssemblage
)
    return Makie.convert_arguments(T, sitestats(asm)[!, stat], asm)
end

function Makie.plot!(p::SiteMap)
    cmap = (
        colormap=p.colormap,
        colorscale=p.colorscale,
        colorrange=p.colorrange,
        lowclip=p.lowclip,
        highclip=p.highclip,
        nan_color=p.nan_color,
        alpha=p.alpha,
    )
    if p.locations[] isa EcoBase.AbstractGridded
        map!(p, [:values, :locations], [:xs, :ys, :image]) do v, grd
            # EcoBase's image is (y, x); Makie's heatmap wants (x, y)
            return (
                collect(xrange(grd, EcoBase.CellCentre())),
                collect(yrange(grd, EcoBase.CellCentre())),
                permutedims(EcoBase.convert_to_image(v, grd)),
            )
        end
        # The sites as colours rather than values, so the colour bar ignores them
        map!(p, [:values, :locations, :empty_color], :sites) do v, grd, c
            issite = permutedims(EcoBase.convert_to_image(ones(length(v)), grd))
            site, other = to_color(c), RGBAf(0, 0, 0, 0)
            return [isnan(x) ? other : site for x in issite]
        end
        heatmap!(p, p.xs, p.ys, p.sites; inspectable=false)
        heatmap!(p, p.xs, p.ys, p.image; cmap...)
    else
        map!(p, :locations, :points) do pnt
            cd = coordinates(pnt, EcoBase.XThenY())
            return Point2d.(cd[:, 1], cd[:, 2])
        end
        scatter!(
            p,
            p.points;
            color=p.empty_color,
            markersize=p.markersize,
            marker=p.marker,
            inspectable=false,
        )
        scatter!(
            p, p.points; color=p.values, markersize=p.markersize, marker=p.marker, cmap...
        )
    end
    return p
end

# Equal scales on both axes. `autolimitaspect` sets the limits from the axis's size, the
# size depends on the tick labels' width and the labels on the limits; for some sizes
# this has no fixed point and Makie recurses until a StackOverflowError (a Makie bug:
# `Axis(Figure(size = (100, 200))[1, 1]; autolimitaspect = 1)` with a heatmap). Tick
# label space that only grows (`:max_auto`) always settles.
const _EQUAL_SCALES = (;
    autolimitaspect=1, xticklabelspace=:max_auto, yticklabelspace=:max_auto
)

# Like EcoBase's `aspect_ratio --> :equal, grid --> false`
function Makie.preferred_axis_attributes(::Type{Axis}, ::SiteMap)
    return (; _EQUAL_SCALES..., xgridvisible=false, ygridvisible=false)
end
