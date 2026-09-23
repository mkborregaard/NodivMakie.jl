# CLAUDE.md — NodivMakie

Private Makie plotting companion to Nodiv (`../Nodiv`, dev'd via `[sources]` path).
Nodiv itself keeps its RecipesBase/Plots recipes; this package reimplements them in Makie
with the goal of interactive figures (clickable tree nodes linked to SOS maps) and
complex Makie layouts.

- Makie 0.24 compute-graph recipes: derive data with `map!(p, inputs, outputs)` inside
  `plot!`; number of child plots is fixed at construction.
- Depends on Makie only; the backend (CairoMakie/GLMakie) is the caller's choice.
  Tests use CairoMakie (`test/Project.toml`, run with `env_path = <root>/test/`).
- Node layout comes from `Phylo._findxy` (internal) so node order matches Phylo and
  Nodiv's `plot_gnd`. Phylo compat is pinned accordingly.
- Do not modify Nodiv or SpatialEcology from here; call out when a change needs them.

## Notes from building the maps and explorer

- Real-data checks use the workshop's cached `res_e`/`res_g` only. Parsing the 9852-tip
  tree is slow; serialise parsed inputs to a scratch dir rather than reparsing.
- CSV 1.x reads names as `DataString`; EcoBase 0.2's `asindices` then cannot match them
  against `String` species names (breaks `get_clade`). Convert `species` to `String`
  when building assemblages (see examples/workshop_explorer.jl). Upstream issue in EcoBase.
- GLMakie cannot open a GL context from Claude's shell on this machine (GLFW segfaults
  in `_glfwGetMonitorPosCocoa`, also unsandboxed). Click handling is tested with CairoMakie
  and a stand-in `pickfn`; real clicks must be checked by the user.
- Clade richness in the node panel uses `CladeRichness` (sites x species index), tested
  equal to `richness(get_clade(...))`; a switch is ~15 ms on the 18k-cell data vs ~0.8 s.
