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
