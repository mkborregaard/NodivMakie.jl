# Which species to show as images around a tree. Tips sit at positions 1..ntips along
# the tree's outer edge (the fan's circle, or the dendrogram's tip column), and every
# clade covers a contiguous run of them. The image size sets how many one-slot-wide
# images fit along the edge; the tree is then cut into as many disjoint clades at least
# one slot wide as possible, and each clade is shown by its widest-ranging species.

"""
    _tip_spans(layout) -> (first, last)

For every node of `layout`, the first and last tip position (1..ntips) of its clade. A
clade's tips are always a contiguous run, as tips are laid out in tree order.
"""
function _tip_spans(l::TreeLayout)
    n = length(l)
    lo, hi = fill(typemax(Int), n), fill(typemin(Int), n)
    for i in 1:n
        l.isleaf[i] || continue
        t = round(Int, l.depth[i])
        j = i
        while j != 0 && (t < lo[j] || t > hi[j])
            lo[j] = min(lo[j], t)
            hi[j] = max(hi[j], t)
            j = l.parent[j]
        end
    end
    return lo, hi
end

"""
    select_clades(layout, nslots; minclade=0.5, circular=false) -> Vector{Int}

Choose the clades (as `layout` indices, in tip order) to show one image each, with the
tip edge divided into `nslots` equal slots of `w = ntips / nslots` tips. An image is one
slot wide and sits at the centre of its clade.

The chosen clades
- are disjoint: no species is in two images;
- are each at least `minclade * w` tips wide, so an image stands for a real part of the
  tree rather than a few species;
- have centres at least `w` tips apart, so neighbouring images never overlap (with
  `circular = true`, as for a fan, also across the gap where the ring closes).

Among the choices that meet these, the one with the most images is taken, and among
those the one covering the most species. So a clade too small to have an image goes with
its sister into their parent's image where that costs no image, and is otherwise left
as a gap. With `minclade = 1` there can be no more images than slots, and every clade
fills at least its own slot.

The choice is exact: clades are taken in tip order, and each is joined to the best
choice among the clades that end before it and whose centre is at least `w` before its
own (a weighted interval scheduling over all clades).
"""
function select_clades(l::TreeLayout, nslots::Integer; minclade=0.5, circular=false)
    lo, hi = _tip_spans(l)
    T = _ntips(l)
    w = T / clamp(nslots, 1, T)
    cand = findall(i -> hi[i] - lo[i] + 1 >= minclade * w - 1e-9, eachindex(lo))
    isempty(cand) && return Int[]
    # Prefix maxima of (images, species covered, clade) over twice the clade centres, in a
    # Fenwick tree
    K = 2T
    fenwick = fill((0, 0, 0), K)
    function raise!(k, v)
        while k <= K
            fenwick[k] = max(fenwick[k], v)
            k += k & -k
        end
        return nothing
    end
    function prefixmax(k)
        v = (0, 0, 0)
        k = min(k, K)
        while k > 0
            v = max(v, fenwick[k])
            k -= k & -k
        end
        return v
    end
    best = Dict{Int,Tuple{Int,Int}}()
    prev = Dict{Int,Int}()
    bylo, byhi = sort(cand; by=i -> lo[i]), sort(cand; by=i -> hi[i])
    next = 1
    for i in bylo
        while next <= length(byhi) && hi[byhi[next]] < lo[i]
            j = byhi[next]
            raise!(lo[j] + hi[j], (best[j]..., j))
            next += 1
        end
        n, cov, j = prefixmax(floor(Int, lo[i] + hi[i] - 2w + 1e-9))
        best[i] = (n + 1, cov + hi[i] - lo[i] + 1)
        prev[i] = j
    end
    c = argmax(i -> (best[i]..., -i), cand)
    chosen = Int[]
    while c != 0
        push!(chosen, c)
        c = prev[c]
    end
    reverse!(chosen)
    # a fan's ring closes: the first and last images must also be a slot apart
    centre(i) = (lo[i] + hi[i]) / 2
    while circular &&
          length(chosen) >= 2 &&
          (T + 1) - (centre(chosen[end]) - centre(chosen[1])) < w - 1e-9
        # drop the smaller of the two
        firstsmaller = hi[chosen[1]] - lo[chosen[1]] <= hi[chosen[end]] - lo[chosen[end]]
        deleteat!(chosen, firstsmaller ? 1 : length(chosen))
    end
    return chosen
end

"""
    CladeImage

One image position around a tree.
- `clade`: the node whose clade the image stands for
- `tips`: the tip positions of the clade
- `species`: the clade's species with the largest range
- `shown`: the species whose image is shown: the largest-range species that has an
  image, or `nothing` if none of the clade's species has one
"""
struct CladeImage
    clade::String
    tips::UnitRange{Int}
    species::String
    shown::Union{String,Nothing}
end

# Range size per species name, from a Dict or an assemblage's occupancy (the number of
# cells each species occupies)
_range_sizes(d::AbstractDict) = Dict(String(k) => Float64(v) for (k, v) in d)
function _range_sizes(asm::EcoBase.AbstractAssemblage)
    return Dict(
        String(n) => Float64(o) for (n, o) in zip(EcoBase.thingnames(asm), occupancy(asm))
    )
end

# Species by range size, largest first; ties alphabetically
_by_range(species, rs) = sort(species; by=s -> (-get(rs, s, 0.0), s))

"""
    image_clades(layout, nslots, rangesize, images; minclade, circular) -> Vector{CladeImage}

The clades chosen by [`select_clades`](@ref) (which takes the keyword arguments), each
with its representative species: the one with the largest range size among those with an
image in `images` (anything with `haskey`, e.g. a [`SpeciesImages`](@ref)). `rangesize` is
a Dict of species => range size or an assemblage (range size = number of occupied cells).
Ties go alphabetically.
"""
function image_clades(l::TreeLayout, nslots::Integer, rangesize, images; kwargs...)
    rs = _range_sizes(rangesize)
    lo, hi = _tip_spans(l)
    tips = findall(l.isleaf)
    tipnames = l.names[tips[sortperm(l.depth[tips])]]    # by position along the tip edge
    return map(select_clades(l, nslots; kwargs...)) do i
        sps = _by_range(tipnames[lo[i]:hi[i]], rs)
        k = findfirst(s -> haskey(images, s), sps)
        return CladeImage(
            l.names[i], lo[i]:hi[i], first(sps), k === nothing ? nothing : sps[k]
        )
    end
end

"""
    image_geometry(layout, treetype; imagesize=automatic, nimages=automatic,
                  gap = 0.04, spacing = 0.1, shape = :circle)

How large the images are and how many fit, in the tree's data coordinates.

For a `:fan`, images sit in a ring just outside the tips. `imagesize` is the image
width as a fraction of the tree's radius (default 0.15), and the ring's centre radius
is `radius * (1 + gap) + width / 2`. The number of images is the most that fit around
the ring with `spacing` (a fraction of the width) between neighbours. Square images
need √2 more room than discs, as upright squares on a slant would otherwise touch.

For a `:dendrogram`, images stand in a column right of the tips. `imagesize` is the
image height as a fraction of the tip column (default 0.05, in tip units), and they
are stacked with `spacing` between them.

Passing `nimages` instead fixes the number, and the size follows from it.

Returns `(; nimages, size, radius)`: `size` in data units (tip units for a dendrogram),
`radius` the fan ring's centre radius (`NaN` for a dendrogram).
"""
function image_geometry(
    l::TreeLayout,
    treetype::Symbol;
    imagesize=automatic,
    nimages=automatic,
    gap=0.04,
    spacing=0.1,
    shape=:circle,
)
    T = _ntips(l)
    room = (1 + spacing) * (shape === :circle ? 1.0 : sqrt(2))
    if treetype === :fan
        R = maximum(l.height)
        span = 2pi * T / (T + 1)                     # angle covered by the tips
        if nimages === automatic
            s = (imagesize === automatic ? 0.15 : imagesize) * R
            r = R * (1 + gap) + s / 2
            n = floor(Int, span / (2asin(min(1, s * room / 2r))))
        else
            # size for n images: the chord between neighbours is room * s at
            # r = R(1 + gap) + s/2; solve by bisection (the chord grows with s)
            n = nimages
            Δ = span / n
            fits(s) = 2 * (R * (1 + gap) + s / 2) * sin(Δ / 2) >= room * s
            a, b = 0.0, 2R
            for _ in 1:60
                m = (a + b) / 2
                if fits(m)
                    a = m
                else
                    b = m
                end
            end
            s = a
            r = R * (1 + gap) + s / 2
        end
        return (; nimages=clamp(n, 1, T), size=s, radius=r)
    elseif treetype === :dendrogram
        if nimages === automatic
            s = (imagesize === automatic ? 0.05 : imagesize) * T
            n = floor(Int, T / (s * room))
        else
            n = nimages
            s = T / (n * room)
        end
        return (; nimages=clamp(n, 1, T), size=s, radius=NaN)
    else
        throw(_treetype_error(treetype))
    end
end
