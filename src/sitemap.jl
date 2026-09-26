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
Missing and NaN values are drawn with `nan_color`, transparent by default.
"""
@recipe SiteMap (values, locations) begin
    "Marker size for point sites."
    markersize = 6
    "Marker shape for point sites."
    marker = :rect
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
    T::Type{<:SiteMap}, stat::Symbol, asm::SpatialEcology.SEAssemblage
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
        heatmap!(p, p.xs, p.ys, p.image; cmap...)
    else
        map!(p, :locations, :points) do pnt
            cd = coordinates(pnt, EcoBase.XThenY())
            return Point2d.(cd[:, 1], cd[:, 2])
        end
        scatter!(
            p, p.points; color=p.values, markersize=p.markersize, marker=p.marker, cmap...
        )
    end
    return p
end

# Like EcoBase's `aspect_ratio --> :equal, grid --> false`
function Makie.preferred_axis_attributes(::Type{Axis}, ::SiteMap)
    return (; autolimitaspect=1, xgridvisible=false, ygridvisible=false)
end
