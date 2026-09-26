# The tree, assemblage and analysis result the tests share

const NEWICK = "((a:1,b:1)n1:1,(c:1.5,(d:1,e:1)n3:0.5)n2:1)root;"

child(p, T) = filter(c -> c isa T, p.plots)
branchlines(p) = only(child(p, Lines))
markers(p) = child(p, Scatter)[2:end]   # the first scatter is the invisible padding

tree = parsenewick(NEWICK)

# Five species on a 4 x 3 grid; `sitecolumns` = one column per site
grid = [(x, y) for y in 1:3 for x in 1:4]
occ = zeros(Int, 5, 12)
occ[1, 1:6] .= 1; occ[2, 4:9] .= 1; occ[3, 7:12] .= 1; occ[4, [1, 5, 9]] .= 1; occ[5, 10:12] .= 1
sites = ["s$i" for i in 1:12]
asm = Assemblage(occ, Float64[first.(grid) last.(grid)], sites, ["a", "b", "c", "d", "e"])
internal = ["root", "n1", "n2", "n3"]                # the internal nodes in tree order
sos = Dict(n => collect(range(-8, 8; length = 12)) .* k for (k, n) in enumerate(internal))
res = NodeAnalysis(internal, Dict("root" => 0.2, "n1" => 0.4, "n2" => 0.9, "n3" => 0.5), sos)
