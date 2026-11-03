module IntegralResurgence

    using DynamicPolynomials, MultivariatePolynomials, LinearAlgebra, GenericSchur, HomotopyContinuation,
        SpecialFunctions, QuadGK
    import Singular_jll, Polynomials

    include("Saddles.jl")
    include("GaussManinSystem.jl")
    include("Asymptotic.jl")
    include("Pade.jl")
    include("Borel-Pade.jl")
    include("Plots.jl")
    include("Hyperlogarithm.jl")
    include("Mould.jl")
    include("SouthEast.jl")
    include("Superasymptotics.jl")
    include("Hyperterminants.jl")

    # Saddle points and the Gauss-Manin system
    export findAllSaddles, refineSaddle, GaussManin

    # Adjacency: Stokes constants from the Gauss-Manin connection
    export stokes_constant, stokes_constants

    # Adjacency from Borel-Padé (numerical alternative)
    export Expansion, asymptotic, findPoles, findStokesConstants, BorelPadeResum

    # Relevance: South-East rule
    export southeast, relevance

    # Asymptotics of the integral
    export transseries, singulants, divergence_rate, superasymptotic
    export hyperterminant, hyperasymptotic, path_lengths

    # Plots
    export plotStokesGraphGM, plotAdjacencyGM, plotAdjacencyPoles

    # Calling a plotting function without a Makie backend gives a MethodError; point to the extension.
    function __init__()
        Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
            if exc.f isa Function && nameof(exc.f) in MAKIE_FUNCTIONS && isnothing(Base.get_extension(IntegralResurgence, :IntegralResurgenceMakieExt))
                print(io, "\n\n", nameof(exc.f), " needs a Makie backend. Load one first, e.g. `using CairoMakie`.")
            end
        end
    end

end
