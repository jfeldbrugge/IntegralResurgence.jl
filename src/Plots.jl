# Plotting. The methods live in the package extension ext/IntegralResurgenceMakieExt.jl, which loads
# with Makie: run `using CairoMakie` (or GLMakie) to enable them. Without a Makie backend the functions
# exist but have no methods, and calling them prints a hint (see `__init__`).

"""
    plotPoles(poles, transseries, plotRange = 3)

Plot the Padé poles of the Borel transforms (`findPoles`), shifted to the Borel singularity of
each series, together with the singularities themselves. Returns `(fig, ax)`.
Requires a Makie backend: load `CairoMakie` (or `GLMakie`) first.
"""
function plotPoles end

"""
    plotAdjacencyPoles(transseries, K, poles, plotRange = 3)

Plot the adjacency graph (an edge for every K[i, j] ≠ 0) of the expansions `transseries`
together with their shifted Padé poles (`findPoles`). Returns `(fig, ax)`.
Requires a Makie backend: load `CairoMakie` (or `GLMakie`) first.
"""
function plotAdjacencyPoles end

"""
    plotAdjacencyGM(Σ, K, plotRange = 3)

Plot the adjacency graph of the Borel singularities Σ: an edge for every K[i, j] ≠ 0, without
the values of the Stokes constants (see `plotStokesGraphGM`). Returns `(fig, ax)`.
Requires a Makie backend: load `CairoMakie` (or `GLMakie`) first.
"""
function plotAdjacencyGM end

"""
    plotStokesGraphGM(Σ, K; plotRange = nothing, from = :column, orient = :positive, atol = 1e-6, ...)

Plot the Stokes graph in the Borel plane: the singularities Σ as points and one edge for every
adjacent pair, labelled by its Stokes constant. K is antisymmetric (K_{i←j} = -K_{j←i}), so an
edge carries a single number once it is given a direction.

- `from = :column` reads K[j, i] as K_{j←i}, the constant of the edge i → j (convention of
  `stokes_constants`); `from = :row` reads K[i, j] as the edge i → j (as in `findStokesConstants`).
- `orient = :positive`: each edge is drawn as an arrow i → j along the direction in which the
  Stokes constant is positive, labelled by |K|. The sign is carried by the arrow.
- `orient = :index`: each edge points from the lower to the higher index, i → j with i < j. It is
  labelled by the signed constant, and the colour gives the sign (`poscolor`, `negcolor`;
  `cplxcolor` when the constant is not real).
- The line width grows with |K|. A warning is issued if K is not antisymmetric within `atol`.
- Each singularity is labelled in the widest angular gap between its edges, so the labels do
  not overlap the graph.
- `plotRange = nothing` fits the limits to Σ; a number r gives the square [-r, r]².

Returns `(fig, ax)`.
Requires a Makie backend: load `CairoMakie` (or `GLMakie`) first.
"""
function plotStokesGraphGM end

"""
Format a Stokes constant: integers without decimals, otherwise rounded to `digits`.
"""
function stokes_label(k; digits = 2, atol = 1e-6)
    if abs(imag(k)) < atol
        r = real(k)
        return abs(r - round(r)) < atol ? string(round(Int, r)) : string(round(r; digits = digits))
    end
    return string(round(ComplexF64(k); digits = digits))
end

"""
Direction (angle) for the label of point i: the bisector of the widest angular gap between the
edges to `nbrs` (or, without edges, between the directions to all other points).
"""
function label_angle(z, i, nbrs)
    others = isempty(nbrs) ? [k for k in eachindex(z) if k != i] : nbrs
    isempty(others) && return 3π / 4
    φ = sort([mod(angle(z[k] - z[i]), 2π) for k in others])
    gaps = [diff(φ); φ[1] + 2π - φ[end]]
    m = argmax(gaps)
    return φ[m] + gaps[m] / 2
end

"Names of the plotting functions provided by the Makie extension."
const MAKIE_FUNCTIONS = (:plotPoles, :plotAdjacencyPoles, :plotAdjacencyGM, :plotStokesGraphGM)
