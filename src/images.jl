# Species images supplied by the user: a directory of files named by species. Nothing
# is bundled with the package, as images are often copyrighted.

"""
    speciesname(s)

The key species names and image file names are matched on: lower case, with runs of
spaces, hyphens and underscores as one underscore. So `Carduelis_hornemanni`,
`Carduelis hornemanni` and `carduelis-hornemanni` all match.
"""
speciesname(s) = lowercase(replace(strip(String(s)), r"[\s_\-]+" => "_"))

"""
    SpeciesImages(dir; extensions = (".png", ".jpg", ".jpeg"), maxpixels = 256)

The species images in `dir`, one file per species, named by species (see
[`speciesname`](@ref) for the matching). Images are loaded when first used, reduced to
at most `maxpixels` on a side, and cached.

`haskey(images, species)` tells whether a species has an image, and `images[species]`
gives it as a matrix of colours. Links whose target is missing are ignored.
"""
struct SpeciesImages
    dir::String
    files::Dict{String, String}          # speciesname => file
    maxpixels::Int
    cache::Dict{String, Matrix{RGBAf}}
end

function SpeciesImages(dir::AbstractString; extensions = (".png", ".jpg", ".jpeg"),
                       maxpixels = 256)
    isdir(dir) || throw(ArgumentError("No such directory: $dir"))
    files = Dict{String, String}()
    exts = lowercase.(collect(extensions))
    # sorted, so a species with several files always gets the same one
    for f in sort(readdir(dir))
        stem, ext = splitext(f)
        lowercase(ext) in exts || continue
        isfile(joinpath(dir, f)) || continue        # e.g. a link to a file not there (yet)
        get!(files, speciesname(stem), joinpath(dir, f))
    end
    return SpeciesImages(String(dir), files, maxpixels, Dict{String, Matrix{RGBAf}}())
end

Base.haskey(imgs::SpeciesImages, species) = haskey(imgs.files, speciesname(species))
Base.length(imgs::SpeciesImages) = length(imgs.files)
Base.show(io::IO, imgs::SpeciesImages) =
    print(io, "SpeciesImages(", length(imgs), " images in ", repr(imgs.dir), ")")

function Base.getindex(imgs::SpeciesImages, species)
    key = speciesname(species)
    haskey(imgs.files, key) || throw(KeyError(species))
    return get!(imgs.cache, key) do
        thumbnail(RGBAf.(FileIO.load(imgs.files[key])), imgs.maxpixels)
    end
end

# Reduce an image to at most `maxpixels` on its longest side by averaging blocks of
# pixels; avoids drawing (and saving) full-size photographs.
function thumbnail(img::AbstractMatrix{RGBAf}, maxpixels)
    f = cld(maximum(size(img)), maxpixels)
    f <= 1 && return Matrix(img)
    h, w = fld.(size(img), f)
    return [sum(@view img[(i - 1) * f + 1:i * f, (j - 1) * f + 1:j * f]) / f^2
            for i in 1:h, j in 1:w]
end

# The central square of an image
function squarecrop(img::AbstractMatrix)
    h, w = size(img)
    n = min(h, w)
    i0, j0 = (h - n) ÷ 2, (w - n) ÷ 2
    return img[i0 + 1:i0 + n, j0 + 1:j0 + n]
end

# Pixels outside the inscribed circle made transparent, with a one-pixel soft edge
function discmask(img::AbstractMatrix{RGBAf})
    n = size(img, 1)
    c = (n + 1) / 2
    r = n / 2
    return [begin
                d = hypot(i - c, j - c)
                a = clamp(r - d + 0.5, 0, 1)
                p = img[i, j]
                RGBAf(p.r, p.g, p.b, p.alpha * a)
            end for i in 1:n, j in 1:size(img, 2)]
end

# The image centred on a transparent n x n canvas, so nothing is cut off
function padto(img::AbstractMatrix{RGBAf}, n)
    h, w = size(img)
    out = fill(RGBAf(0, 0, 0, 0), n, n)
    i0, j0 = (n - h) ÷ 2, (n - w) ÷ 2
    out[i0 + 1:i0 + h, j0 + 1:j0 + w] .= img
    return out
end

# Near-white pixels made transparent (with a soft edge), for illustrations on white
function keywhite(img::AbstractMatrix{RGBAf}; threshold = 0.92)
    return map(img) do p
        w = min(p.r, p.g, p.b)                      # how white: the darkest channel
        a = clamp((1 - w) / (1 - threshold), 0, 1)
        RGBAf(p.r, p.g, p.b, p.alpha * a)
    end
end

# The image centred on a transparent canvas just large enough for its visible content:
# for `shape = :circle` the canvas is the circle around the content, so the content fills
# the disc. The content is the pixels at least half opaque (faint ones, e.g. off-white
# specks left after keying out a white background, are drawn but do not set the size);
# for an image without transparency it is the whole image.
#
# `clip` is the share of the content allowed outside the disc: 0 keeps all of it (the
# circle is centred on the content's bounding box or its centre of mass, whichever needs
# the smaller circle); e.g. 0.02 lets the outermost 2% (a tail tip, the feet) be cut off
# for a larger image, with the circle centred on the content's centre of mass.
function fitcontent(img::AbstractMatrix{RGBAf}, shape::Symbol; margin = 0.03, clip = 0.0)
    idx = findall(p -> p.alpha >= 0.5, img)
    isempty(idx) && return padto(img, maximum(size(img)))
    h, w = size(img)
    if shape === :circle
        (r1, r2), (c1, c2) = extrema(i[1] for i in idx), extrema(i[2] for i in idx)
        centres = [((r1 + r2) / 2, (c1 + c2) / 2),
                   (sum(i[1] for i in idx) / length(idx), sum(i[2] for i in idx) / length(idx))]
        clip > 0 && (centres = centres[2:2])
        best = (Inf, 0.0, 0.0)
        for (ci, cj) in centres
            d = sort!([hypot(i[1] - ci, i[2] - cj) for i in idx])
            reach = d[clamp(ceil(Int, (1 - clip) * length(d)), 1, length(d))] + 1
            reach < best[1] && (best = (reach, ci, cj))
        end
        reach, ci, cj = best
        n = ceil(Int, 2reach * (1 + margin))
        # place the image so the centre lands mid-canvas; parts outside the canvas are cut
        out = fill(RGBAf(0, 0, 0, 0), n, n)
        oi, oj = round(Int, (n + 1) / 2 - ci), round(Int, (n + 1) / 2 - cj)
        for i in max(1, 1 - oi):min(h, n - oi), j in max(1, 1 - oj):min(w, n - oj)
            out[i + oi, j + oj] = img[i, j]
        end
        return discmask(out)
    end
    (r1, r2), (c1, c2) = extrema(i[1] for i in idx), extrema(i[2] for i in idx)
    sub = img[r1:r2, c1:c2]
    return padto(sub, ceil(Int, maximum(size(sub)) * (1 + margin)))
end

# An image ready to draw with `image!` over a square, turned so it is upright (`image!`
# puts the first matrix index along x and the second along y, upwards):
# - `fit = :pad` shows the whole image, as large as its circle or square allows; with a
#   transparent background only the visible content has to fit, so it fills the disc
# - `fit = :crop` fills the circle or square with the image's central part
# - `whitebackground = true` makes a white background transparent first
# - `clip`: with `:pad` on a disc, the share of the content allowed to be cut off at the
#   disc's edge for a larger image (see `fitcontent`)
function markerimage(img::AbstractMatrix{RGBAf}, shape::Symbol; fit = :pad,
                     whitebackground = false, clip = 0.0)
    whitebackground && (img = keywhite(img))
    if fit === :crop
        sq = squarecrop(img)
        shape === :circle && (sq = discmask(sq))
    else
        sq = fitcontent(img, shape; clip)
    end
    return rotr90(sq)
end
