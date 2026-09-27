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
slot wide and sits over its clade: its centre is at one of the clade's tip positions
(see [`image_positions`](@ref) for where).

The chosen clades
- are disjoint: no species is in two images;
- are each at least `minclade * w` tips wide, so an image stands for a real part of the
  tree rather than a few species;
- can have their images at least `w` tips apart, so neighbouring images never overlap.
  With `circular = true`, as for a fan, images also keep half a slot from the gap where
  the ring closes, so the first and last do not overlap either.

Among the choices that meet these, the one with the most images is taken, and then, as
far as that costs no image, the one covering more species. So a clade too small to have
an image goes with its sister into their parent's image where that costs no image, and
is otherwise left as a gap. With `minclade = 1` there can be no more images than slots,
and every clade fills at least its own slot.

The number of images is exact: clades are taken in tip order, and for each clade and
number of images the earliest position its image can have is kept, from the choices among
the clades that end before it.
"""
function select_clades(l::TreeLayout, nslots::Integer; minclade=0.5, circular=false)
    lo, hi, w = _image_bounds(l, nslots)
    T = _ntips(l)
    cand = findall(eachindex(lo)) do i
        return hi[i] - lo[i] + 1 >= minclade * w - 1e-9 &&
               _lower(lo, i, w, circular) <= _upper(hi, i, T, w, circular)
    end
    isempty(cand) && return Int[]
    K = if circular
        floor(Int, (T + 1 - w) / w + 1e-9) + 1
    else
        floor(Int, (T - 1) / w + 1e-9) + 1
    end
    K = max(K, 1)
    # For clade i with k images ending in it: (earliest position of its image, species
    # covered, previous clade). `fronts[k]` holds the choices of k images among the
    # clades that end before the current one, by position, keeping only those that cover
    # more species than every choice with an earlier position.
    state = Dict{Int,Vector{Tuple{Float64,Int,Int}}}()
    fronts = [Tuple{Float64,Int,Int}[] for _ in 1:K]
    bylo, byhi = sort(cand; by=i -> lo[i]), sort(cand; by=i -> hi[i])
    next = 1
    for i in bylo
        while next <= length(byhi) && hi[byhi[next]] < lo[i]
            j = byhi[next]
            for (k, (p, cov, _)) in enumerate(state[j])
                isfinite(p) && _insert_front!(fronts[k], (p, cov, j))
            end
            next += 1
        end
        a, b = _lower(lo, i, w, circular), _upper(hi, i, T, w, circular)
        n = hi[i] - lo[i] + 1
        s = fill((Inf, 0, 0), K)
        s[1] = (a, n, 0)
        for k in 2:K
            f = fronts[k - 1]
            isempty(f) && break
            # any earlier image at most a - w lets this one be at a: take the one
            # covering the most species; otherwise the earliest
            m = searchsortedlast(f, (a - w + 1e-9, typemax(Int), typemax(Int)))
            p, cov, j = m > 0 ? (a, f[m][2], f[m][3]) : (f[1][1] + w, f[1][2], f[1][3])
            p <= b + 1e-9 && (s[k] = (p, cov + n, j))
        end
        state[i] = s
    end
    kmax = maximum(i -> findlast(e -> isfinite(e[1]), state[i]), cand)
    c = argmax(i -> (isfinite(state[i][kmax][1]) ? state[i][kmax][2] : -1, -i), cand)
    chosen = Int[]
    for k in kmax:-1:1
        push!(chosen, c)
        c = state[c][k][3]
    end
    return reverse!(chosen)
end

# Add (position, species covered, clade) to a front sorted by position, in which later
# positions cover more species
function _insert_front!(f, e)
    m = searchsortedlast(f, e)
    m > 0 && f[m][2] >= e[2] && return f           # an earlier choice covers as many
    insert!(f, m + 1, e)
    while m + 2 <= length(f) && f[m + 2][2] <= e[2]
        deleteat!(f, m + 2)
    end
    return f
end

# Tip spans of all clades and the slot width. The tip positions an image centre can take
# for clade i are _lower(...) .. _upper(...): the clade's tips, and for a ring at least
# half a slot from the gap where it closes (at 0 ≡ ntips + 1).
function _image_bounds(l::TreeLayout, nslots)
    lo, hi = _tip_spans(l)
    T = _ntips(l)
    return lo, hi, T / clamp(nslots, 1, T)
end
_lower(lo, i, w, circular) = circular ? max(Float64(lo[i]), w / 2) : Float64(lo[i])
function _upper(hi, i, T, w, circular)
    return circular ? min(Float64(hi[i]), T + 1 - w / 2) : Float64(hi[i])
end

"""
    image_positions(layout, nslots, clades; circular=false) -> Vector{Float64}

The tip positions of the images for `clades` chosen by [`select_clades`](@ref) (with the
same `nslots` and `circular`). Each image is as near the centre of its clade as the
images beside it allow, and always over one of the clade's tips; images are at least a
slot (`ntips / nslots` tips) apart. Where they are crowded, the images are placed so that
the sum of squared distances to their clades' centres is smallest.
"""
function image_positions(l::TreeLayout, nslots::Integer, clades; circular=false)
    lo, hi, w = _image_bounds(l, nslots)
    T, n = _ntips(l), length(clades)
    n == 0 && return Float64[]
    a = [_lower(lo, i, w, circular) for i in clades]
    b = [_upper(hi, i, T, w, circular) for i in clades]
    # the earliest and latest each image can be
    L, R = copy(a), copy(b)
    for k in 2:n
        L[k] = max(L[k], L[k - 1] + w)
    end
    for k in (n - 1):-1:1
        R[k] = min(R[k], R[k + 1] - w)
    end
    all(k -> L[k] <= R[k] + 1e-9, 1:n) ||
        throw(ArgumentError("the images of these clades cannot be a slot apart"))
    # As near the centres as possible (least squares): with q[k] = p[k] - (k - 1)w the
    # images a slot apart is q nondecreasing, within bounds that are nondecreasing too
    shift = (0:(n - 1)) .* w
    q = _isotonic([(lo[i] + hi[i]) / 2 for i in clades] .- shift, L .- shift, R .- shift)
    return q .+ shift
end

# The nondecreasing x with lower[k] <= x[k] <= upper[k] nearest to y (least squares),
# for nondecreasing bounds that admit one, by pooling adjacent violators: runs of equal
# x form blocks, each at the mean of its y clamped to the bounds of all its members
function _isotonic(y, lower, upper)
    blocks = Tuple{Float64,Int,Float64,Int}[]            # (value, size, sum of y, first)
    for k in eachindex(y)
        push!(blocks, (clamp(y[k], lower[k], upper[k]), 1, y[k], k))
        while length(blocks) > 1 && blocks[end - 1][1] > blocks[end][1]
            _, n2, s2, _ = pop!(blocks)
            _, n1, s1, f = blocks[end]
            # the bounds all members share: nondecreasing, so the last's lower bound and
            # the first's upper bound
            v = clamp((s1 + s2) / (n1 + n2), lower[k], upper[f])
            blocks[end] = (v, n1 + n2, s1 + s2, f)
        end
    end
    return reduce(vcat, [fill(v, n) for (v, n, _, _) in blocks]; init=Float64[])
end

"""
    CladeImage

One image position around a tree.
- `clade`: the node whose clade the image stands for
- `tips`: the tip positions of the clade
- `position`: the tip position of the image's centre (see [`image_positions`](@ref))
- `species`: the clade's species with the largest range
- `shown`: the species whose image is shown: the largest-range species that has an
  image, or `nothing` if none of the clade's species has one
"""
struct CladeImage
    clade::String
    tips::UnitRange{Int}
    position::Float64
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
with its image position from [`image_positions`](@ref) and its representative species:
the one with the largest range size among those with an image in `images` (anything with
`haskey`, e.g. a [`SpeciesImages`](@ref)). `rangesize` is a Dict of species => range size
or an assemblage (range size = number of occupied cells). Ties go alphabetically.
"""
function image_clades(
    l::TreeLayout, nslots::Integer, rangesize, images; minclade=0.5, circular=false
)
    rs = _range_sizes(rangesize)
    lo, hi = _tip_spans(l)
    tips = findall(l.isleaf)
    tipnames = l.names[tips[sortperm(l.depth[tips])]]    # by position along the tip edge
    chosen = select_clades(l, nslots; minclade, circular)
    positions = image_positions(l, nslots, chosen; circular)
    return map(chosen, positions) do i, p
        sps = _by_range(tipnames[lo[i]:hi[i]], rs)
        k = findfirst(s -> haskey(images, s), sps)
        return CladeImage(
            l.names[i], lo[i]:hi[i], p, first(sps), k === nothing ? nothing : sps[k]
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
