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

# An image ready to draw with `image!` over a square, turned so it is upright (`image!`
# puts the first matrix index along x and the second along y, upwards):
# - `fit = :pad` shrinks it to fit the square, whole; `:crop` fills the square with its
#   central part (for a disc, :pad also shrinks it to fit inside the circle)
# - `whitebackground = true` makes a white background transparent
function markerimage(img::AbstractMatrix{RGBAf}, shape::Symbol; fit = :pad,
                     whitebackground = false)
    whitebackground && (img = keywhite(img))
    if fit === :crop
        sq = squarecrop(img)
        shape === :circle && (sq = discmask(sq))
    else
        # whole image; inside a disc its corners must fit in the circle, so the canvas is
        # as wide as the image's diagonal
        h, w = size(img)
        sq = padto(img, shape === :circle ? ceil(Int, hypot(h, w)) : max(h, w))
        shape === :circle && (sq = discmask(sq))
    end
    return rotr90(sq)
end
