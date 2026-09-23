# Open the interactive node explorer on the NodivWorkshop data.
#
# Reads the workshop's cleaned inputs and its CACHED node analysis
# (data/node_analysis.jld2); nothing is recomputed. Run in an environment with
# NodivMakie, GLMakie, CSV, DataFrames and JLD2:
#
#     include("examples/workshop_explorer.jl")
#
# Parsing the ~10k-tip tree takes a while the first time in a session.

using NodivMakie, GLMakie, CSV, DataFrames, JLD2

const WORKSHOP = expanduser("~/github/NodivWorkshop")
workshop(path) = joinpath(WORKSHOP, path)

# Site ids stay strings after CSV. Species names are converted to String because
# EcoBase 0.2 cannot match Nodiv's String names against CSV 1.x's DataString names.
function readclean(file)
    df = CSV.read(workshop("data/clean/$file"), DataFrame)
    df.site = string.(df.site)
    "species" in names(df) && (df.species = String.(df.species))
    return df
end

tree = parsenewick(read(workshop("data/clean/tree.nwk"), String))
birds_e = Assemblage(readclean("phylocom_e.csv"), readclean("coords_e.csv"))
birds_g = Assemblage(readclean("phylocom_g.csv"), readclean("coords_g.csv"))
res_e, res_g = load(workshop("data/node_analysis.jld2"), "res_e", "res_g")

# as in the workshop script: the divergent nodes by RMS-SOS
metric = :rms
divergent_e = divergent_nodes(res_e; by = metric, threshold = 2)
divergent_g = divergent_nodes(res_g; by = metric, threshold = 2)

fig_g, ex_g = nodeexplorer(birds_g, tree, res_g; metric, nodes = divergent_g)
display(GLMakie.Screen(), fig_g)
fig_e, ex_e = nodeexplorer(birds_e, tree, res_e; metric, nodes = divergent_e)
display(GLMakie.Screen(), fig_e)

# Linking the two: show the environmental panel for whatever node is clicked in
# geographic space (if it has an SOS there too).
on(ex_g.panel.node) do n
    hassos(tree, res_e.sos, n) && (ex_e.panel.node[] = n)
end
