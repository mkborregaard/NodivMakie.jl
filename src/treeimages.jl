# Drawing species images around a tree: a ring outside a fan, a column right of a
# dendrogram.

"""
    TreeImages

The species images drawn by [`treeimages!`](@ref).
- `clades`: one [`CladeImage`](@ref) per image position, in tip order
- `geometry`: `(; nimages, size, radius)` from [`imagegeometry`](@ref)
- `axis`: the axis the images are in (for a dendrogram, a narrow axis beside the tree)
- `plots`: the image plots
- `pixelsize`: `Observable` of the images' width on screen, in pixels (it follows the
  figure's size); e.g. for images elsewhere at the same size, see [`cladeimages!`](@ref)
"""
struct TreeImages
    clades::Vector{CladeImage}
    geometry::NamedTuple
    axis::Axis
    plots::Vector{Any}
    pixelsize::Observable{Float64}
end

"""
    missingimages(treeimages) -> Vector{String}

The species to find images for: for every image position where none of the clade's
species has an image, the clade's species with the largest range.
"""
missingimages(ti::TreeImages) = [c.species for c in ti.clades if c.shown === nothing]

# An outline around an image of width `w` centred at (x, y) in axis units `sx`, `sy`
# per unit of width (as a polyline, so it stays a circle on screen however the axis
# is scaled), NaN-separated for drawing many at once
function outlinepoints!(pts, x, y, w, sx, sy, shape)
    if shape === :circle
        for φ in range(0, 2pi; length = 65)
            push!(pts, Point2d(x + w / 2 * sx * cos(φ), y + w / 2 * sy * sin(φ)))
        end
    else
        for (a, b) in ((-1, -1), (1, -1), (1, 1), (-1, 1), (-1, -1))
            push!(pts, Point2d(x + a * w / 2 * sx, y + b * w / 2 * sy))
        end
    end
    push!(pts, Point2d(NaN, NaN))
end

# Whether the cursor is over a visible pixel of an image plot (at least half opaque), so
# the transparent surroundings of a bird do not count as the image
function overimage(plot)
    scene = Makie.parent_scene(plot)
    pos = plot.space[] === :pixel ? Makie.mouseposition_px(scene) : Makie.mouseposition(scene)
    (x0, x1), (y0, y1) = plot[1][].data, plot[2][].data
    img = Makie._to_array(plot[3][])
    u, v = (pos[1] - x0) / (x1 - x0), (pos[2] - y0) / (y1 - y0)
    (0 <= u < 1 && 0 <= v < 1) || return false
    i = clamp(floor(Int, u * size(img, 1)) + 1, 1, size(img, 1))
    j = clamp(floor(Int, v * size(img, 2)) + 1, 1, size(img, 2))
    return alpha(img[i, j]) >= 0.5
end

# Hovering over an image's visible pixels (with a `DataInspector`) shows `text` (a string
# or an Observable of one) as a plain tooltip; over its transparent parts, whatever is
# underneath is inspected instead. Makie's own image tooltip shows the pixel's colour and
# marks the pixel; this shows only the text. The text is also the plot's `inspector_label`.
function imagehover!(p, text)
    p.inspector_label = (plot, idx, pos) -> string(to_value(text))
    p.inspector_hover = function (inspector, plot, idx)
        overimage(plot) || return false
        pos = Point2f(Makie.mouseposition_px(inspector.root))
        Makie.update_tooltip_alignment!(inspector, pos; text = string(to_value(text)))
        return true
    end
    p.inspectable = true
    return p
end

displayname(sp) = replace(sp, "_" => " ")

"""
    treeimages!(ax, treeplot, images, rangesize; kwargs...) -> TreeImages

Draw species images around a `treeplot` in `ax`: a ring outside a fan tree, or a
column right of a dendrogram (in a narrow axis added beside `ax`, linked in y).

- `images`: a [`SpeciesImages`](@ref), or the path of a directory of images named by
  species (e.g. `Carduelis_hornemanni.jpg`). No images come with the package.
- `rangesize`: an assemblage (range size = number of occupied cells) or a Dict of
  species => range size.

How many images fit follows from their size (see [`imagegeometry`](@ref)); that many
equal slots divide the tips. Then as many disjoint clades as possible are chosen, each at
least `minclade` of a slot wide and with centres at least a slot apart, so images never
overlap (see [`selectclades`](@ref)). Each is shown by an image at its centre: its
largest-range species that has an image. Positions where no species has an image are
left empty; use [`missingimages`](@ref) to list the species to find images for.

Keyword arguments:
- `imagesize`, `nimages`, `gap = 0.04`, `spacing = 0.1`: size and number of images, see
  [`imagegeometry`](@ref)
- `minclade = 0.5`: the smallest clade given an image, as a fraction of a slot. Lower
  gives more images standing for smaller clades; 1 makes every clade fill its slot.
- `shape = :circle`: `:circle` or `:square` images
- `fit = :pad`: the whole image, shrunk to fit its circle or square; `:crop` fills the
  circle or square with the image's central part (better for photos, but may cut off
  parts of an illustration)
- `whitebackground = false`: make a white background transparent, for illustrations
  drawn on white. The image is then sized so the subject itself fills the disc.
- `clip = 0`: the share of the subject that may be cut off at the disc's edge for a larger
  image, e.g. 0.02 (a tail tip or the feet). 0 never cuts anything.
- `strokewidth = 0`, `strokecolor = :gray40`: an outline around each image (none by
  default)
- `showclades = false`: draw a thin line along the tips each image stands for, in
  `cladecolor = :gray50`

For a fan, the axis must keep a 1:1 aspect (`treeplot` and `nodeexplorer` set this when
they create the axis). Turn tip labels off (`showtips = false`), as the images take
their place.

Hovering over a bird (with a `DataInspector`) shows its species and the clade it stands
for; its transparent surroundings do not count.
"""
function treeimages!(ax::Axis, tp::TreePlot, images, rangesize; imagesize = automatic,
                     nimages = automatic, gap = 0.04, spacing = 0.1, shape = :circle,
                     strokecolor = :gray40, strokewidth = 0, showclades = false,
                     cladecolor = :gray50, minclade = 0.5, fit = :pad,
                     whitebackground = false, clip = 0.0)
    shape in (:circle, :square) || throw(ArgumentError("`shape` must be :circle or :square"))
    images isa AbstractString && (images = SpeciesImages(images))
    l, tt = tp.tree_layout[], tp.treetype[]
    geo = imagegeometry(l, tt; imagesize, nimages, gap, spacing, shape)
    clades = imageclades(l, geo.nimages, rangesize, images; minclade,
                         circular = tt === :fan)
    nmissing = count(c -> c.shown === nothing, clades)
    nmissing > 0 &&
        @info "$nmissing of $(length(clades)) image positions have no image for any species of their clade; `missingimages` lists the species to find"
    centre(c) = (first(c.tips) + last(c.tips)) / 2
    hovertext(c) = "$(displayname(c.shown))\nfor $(c.clade) ($(length(c.tips)) species)"
    shown = filter(c -> c.shown !== nothing, clades)
    s = geo.size
    plots = Any[]

    if tt === :fan
        T = ntips(l)
        pts = [polar(geo.radius, fanangle(centre(c), T)) for c in shown]
        for (c, p) in zip(shown, pts)
            img = markerimage(images[c.shown], shape; fit, whitebackground, clip)
            ip = image!(ax, (p[1] - s / 2) .. (p[1] + s / 2), (p[2] - s / 2) .. (p[2] + s / 2), img)
            push!(plots, imagehover!(ip, hovertext(c)))
        end
        outline = Point2d[]
        for p in pts
            outlinepoints!(outline, p[1], p[2], s, 1, 1, shape)
        end
        imgax = ax
        pixelsize = lift((vp, lims) -> s * vp.widths[1] / max(lims.widths[1], eps()),
                         ax.scene.viewport, ax.finallimits)
        if showclades
            R = maximum(l.height)
            arcs = Point2d[]
            for c in clades
                θs = range(fanangle(first(c.tips) - 0.4, T), fanangle(last(c.tips) + 0.4, T);
                           length = 2 + ceil(Int, 60 * length(c.tips) / T))
                append!(arcs, polar.(R * (1 + gap / 2), θs))
                push!(arcs, Point2d(NaN, NaN))
            end
            push!(plots, lines!(ax, arcs; color = cladecolor, linewidth = 1.5, inspectable = false))
        end
    else
        # A narrow axis right of the tree, sharing its y axis. Its x axis runs from
        # -gap to 1 across the image column, and its width is set so that an image one
        # unit wide is as wide on screen as it is tall: images stay square however the
        # tree axis is shaped.
        gc = ax.layoutobservables.gridcontent[]
        gl, rows, col = gc.parent, gc.span.rows, gc.span.cols.stop + 1
        col <= ncols(gl) && insertcols!(gl, col, 1)
        g = 4gap
        imgax = Axis(gl[rows, col]; limits = ((-g, 1.0), nothing), xgridvisible = false,
                     ygridvisible = false)
        hidedecorations!(imgax)
        hidespines!(imgax)
        linkyaxes!(ax, imgax)
        pixels_per_tip = lift((vp, lims) -> vp.widths[2] / max(lims.widths[2], eps()),
                              ax.scene.viewport, ax.finallimits)
        on(pxy -> (imgax.width = max(1.0, s * pxy * (1 + g))), pixels_per_tip; update = true)
        pixelsize = lift(pxy -> s * pxy, pixels_per_tip)
        for c in shown
            img = markerimage(images[c.shown], shape; fit, whitebackground, clip)
            y = centre(c)
            ip = image!(imgax, 0.0 .. 1.0, (y - s / 2) .. (y + s / 2), img)
            push!(plots, imagehover!(ip, hovertext(c)))
        end
        outline = Point2d[]
        for c in shown
            outlinepoints!(outline, 0.5, centre(c), s, 1 / s, 1, shape)
        end
        if showclades
            H = maximum(l.height)
            bars = Point2d[]
            for c in clades
                push!(bars, Point2d(H * 1.02, first(c.tips) - 0.4),
                      Point2d(H * 1.02, last(c.tips) + 0.4), Point2d(NaN, NaN))
            end
            push!(plots, lines!(ax, bars; color = cladecolor, linewidth = 1.5, inspectable = false))
        end
    end
    if strokewidth > 0 && !isempty(outline)
        push!(plots, lines!(imgax, outline; color = strokecolor, linewidth = strokewidth,
                            inspectable = false))
    end
    return TreeImages(clades, geo, imgax, plots, pixelsize)
end

"""
    cladeimages!(panel, tree, images, rangesize; pixelsize = 60, kwargs...) -> Vector

Put a species image in the top-right corner of the two child-clade richness maps of a
[`NodePanel`](@ref): each clade's species with the largest range size that has an image.
The images follow the node shown, and a clade with no image gets none.

- `images`, `rangesize`: as for [`treeimages!`](@ref)
- `pixelsize = 60`: the images' width on screen, in pixels, or an `Observable` of it (the
  explorer passes its tree images' `pixelsize`, so both are the same size)
- `margin = 6`: the gap to the map's corner, in pixels
- `shape`, `fit`, `whitebackground`, `clip`: as for [`treeimages!`](@ref)

The images are drawn in the maps' screen space: they keep their place and size when the
maps are zoomed or panned. Hovering over the bird (with a `DataInspector`) shows its
species.
Returns the two image plots.
"""
function cladeimages!(np::NodePanel, tree, images, rangesize; pixelsize = 60, margin = 6,
                      shape = :circle, fit = :pad, whitebackground = false, clip = 0.0)
    images isa AbstractString && (images = SpeciesImages(images))
    rs = rangesizes(rangesize)
    function representative(clade)
        sps = filter(s -> haskey(images, s), nodespecies(tree, clade))
        isempty(sps) && return nothing
        return first(sort(sps; by = s -> (-get(rs, s, 0.0), s)))
    end
    child(n, k) = getnodename(tree, getchildren(tree, n)[k])
    px = pixelsize isa Observable ? pixelsize : Observable(Float64(pixelsize))
    blank = fill(RGBAf(0, 0, 0, 0), 2, 2)
    plots = Any[]
    for k in 1:2
        ax = np.axes[2 + k]
        species = Observable{Union{String, Nothing}}(representative(child(np.node[], k)))
        on(n -> (species[] = representative(child(n, k))), np.node)
        img = lift(sp -> sp === nothing ? blank :
                         markerimage(images[sp], shape; fit, whitebackground, clip), species)
        xs = lift((vp, s) -> (vp.widths[1] - margin - s) .. (vp.widths[1] - margin),
                  ax.scene.viewport, px)
        ys = lift((vp, s) -> (vp.widths[2] - margin - s) .. (vp.widths[2] - margin),
                  ax.scene.viewport, px)
        p = image!(ax.scene, xs, ys, img; space = :pixel, visible = lift(!isnothing, species))
        imagehover!(p, lift(sp -> sp === nothing ? "" : displayname(sp), species))
        translate!(p, 0, 0, 1)                   # above the map
        push!(plots, p)
    end
    return plots
end
