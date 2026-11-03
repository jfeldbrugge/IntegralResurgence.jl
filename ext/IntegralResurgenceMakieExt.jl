module IntegralResurgenceMakieExt

# Makie methods of the plotting functions declared in IntegralResurgence (src/Plots.jl). The extension
# loads with Makie, i.e. as soon as the user loads a backend such as CairoMakie or GLMakie.

using IntegralResurgence, Makie
import IntegralResurgence: plotPoles, plotAdjacencyPoles, plotAdjacencyGM, plotStokesGraphGM
using IntegralResurgence: stokes_label, label_angle

function plotPoles(poles, transseries, plotRange = 3)
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "Re", ylabel = "Im", limits=(-plotRange, plotRange,-plotRange, plotRange), aspect = 1)

    for i in eachindex(transseries)
        points = poles[i] .+ transseries[i].f₀
        scatter!(ax, real.(points), imag.(points), label = string(i))
    end

    for i in eachindex(transseries)
        f0 = transseries[i].f₀
        scatter!(ax, [real(f0)], [imag(f0)], color = :black)
        text!(ax, string(i), position = Point2f(real(f0), imag(f0)), align = (:right, :bottom), offset = (8, 8))
    end
    axislegend(ax)
    return fig, ax
end

function plotAdjacencyPoles(transseries, K, poles, plotRange = 3)
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "Re", ylabel = "Im", limits=(-plotRange, plotRange,-plotRange, plotRange), aspect = 1)

    for i in eachindex(transseries)
        points = poles[i] .+ transseries[i].f₀
        scatter!(ax, real.(points), imag.(points), label = string(i))
    end

    for i in axes(K, 1), j in axes(K, 2)
        if K[i, j] != 0
            f0i = transseries[i].f₀
            f0j = transseries[j].f₀

            lines!(ax, [Point2f(real(f0i), imag(f0i)), Point2f(real(f0j), imag(f0j))], color = :red, linewidth = 2)
        end
    end

    for i in eachindex(transseries)
        f0 = transseries[i].f₀
        scatter!(ax, [real(f0)], [imag(f0)], color = :black)
        text!(ax, string(i), position = Point2f(real(f0), imag(f0)), align = (:right, :bottom), offset = (8, 8))
    end

    axislegend(ax)
    return fig, ax
end

function plotAdjacencyGM(Σ, K, plotRange = 3)
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "Re", ylabel = "Im", limits=(-plotRange, plotRange,-plotRange, plotRange), aspect = 1)

    for i in axes(K, 1), j in axes(K, 2)
        if K[i, j] != 0
            f0i = Σ[i]
            f0j = Σ[j]

            lines!(ax, [Point2f(real(f0i), imag(f0i)), Point2f(real(f0j), imag(f0j))], color = :red, linewidth = 2)
        end
    end

    for i in eachindex(Σ)
        f0 = Σ[i]
        scatter!(ax, [real(f0)], [imag(f0)], color = :black)
        text!(ax, string(i), position = Point2f(real(f0), imag(f0)), align = (:right, :bottom), offset = (8, 8), label = string(i))
    end

    axislegend(ax)
    return fig, ax
end

function plotStokesGraphGM(Σ, K; plotRange = nothing, from = :column, orient = :positive, atol = 1e-6,
        edgecolor = :red, poscolor = :royalblue, negcolor = :crimson, cplxcolor = :gray50,
        linewidth = 2, labelsize = 13, labeloffset = 9, arrowsize = 14, digits = 2, legend = true)
    n = length(Σ)
    size(K) == (n, n) || error("plotStokesGraphGM: K must be $(n)×$(n)")
    from in (:column, :row) || error("plotStokesGraphGM: `from` must be :column or :row")
    orient in (:positive, :index) || error("plotStokesGraphGM: `orient` must be :positive or :index")
    Kc = from === :column ? K : transpose(K)          # Kc[j, i] = constant of the edge i → j
    maximum(abs, Kc + transpose(Kc); init = 0.0) < atol ||
        @warn "plotStokesGraphGM: K is not antisymmetric; using the entries with i < j"
    z = ComplexF64.(Σ)

    fig = Figure()
    if plotRange === nothing
        xs, ys = real.(z), imag.(z)
        w = max(maximum(xs) - minimum(xs), maximum(ys) - minimum(ys), 1.0)
        cx, cy = (maximum(xs) + minimum(xs)) / 2, (maximum(ys) + minimum(ys)) / 2
        lims = (cx - 0.65w, cx + 0.65w, cy - 0.65w, cy + 0.65w)
    else
        lims = (-plotRange, plotRange, -plotRange, plotRange)
    end
    ax = Axis(fig[1, 1], xlabel = "Re t", ylabel = "Im t", limits = lims, aspect = DataAspect())

    kmax = maximum(abs, Kc; init = 0.0)
    nbrs = [Int[] for _ in 1:n]                       # neighbours in the drawn graph
    P(z) = Point2f(real(z), imag(z))
    isreal_(k) = abs(imag(k)) < atol

    for i in 1:n, j in i+1:n
        k = Kc[j, i]                                   # K_{j←i}, edge i → j
        abs(k) > atol || continue
        a, b = i, j
        if orient === :positive
            isreal_(k) && real(k) < 0 && ((a, b, k) = (j, i, -k))
            col = isreal_(k) ? edgecolor : cplxcolor
            lab = stokes_label(k; digits = digits, atol = atol)
        else
            col = isreal_(k) ? (real(k) > 0 ? poscolor : negcolor) : cplxcolor
            lab = stokes_label(k; digits = digits, atol = atol)
        end
        push!(nbrs[a], b); push!(nbrs[b], a)
        p0, p1 = z[a], z[b]
        lines!(ax, [P(p0), P(p1)], color = col, linewidth = linewidth * (1 + abs(k) / max(kmax, atol)))
        scatter!(ax, [P(p0 + 0.72 * (p1 - p0))], marker = :rtriangle, markersize = arrowsize,
            rotation = angle(p1 - p0), color = col)
        pm = P((p0 + p1) / 2)                          # label disc at the midpoint
        scatter!(ax, [pm], marker = :circle, markersize = 2.2labelsize, color = :white,
            strokecolor = col, strokewidth = 1.5)
        text!(ax, lab, position = pm, align = (:center, :center), fontsize = labelsize, color = :black)
    end

    scatter!(ax, real.(z), imag.(z), color = :black, markersize = 10)
    for i in 1:n                                       # label in the widest gap between the edges
        θ = label_angle(z, i, nbrs[i])
        c, s = cos(θ), sin(θ)
        text!(ax, string(i), position = P(z[i]), offset = labeloffset .* (c, s),
            align = (c > 0.4 ? :left : c < -0.4 ? :right : :center, s > 0.4 ? :bottom : s < -0.4 ? :top : :center),
            fontsize = labelsize + 2, font = :bold, color = :black, glowcolor = (:white, 0.9), glowwidth = 2)
    end

    if legend && orient === :index
        elems = [LineElement(color = poscolor, linewidth = 3), LineElement(color = negcolor, linewidth = 3)]
        names = Any[L"K_{j \leftarrow i} > 0", L"K_{j \leftarrow i} < 0", L"K \notin \mathbb{R}"]
        if any(k -> abs(k) > atol && !isreal_(k), Kc)
            push!(elems, LineElement(color = cplxcolor, linewidth = 3))
        end
        axislegend(ax, elems, names[1:length(elems)], position = :rb)
    end
    return fig, ax
end

end
