# Hyperlogarithms (src/Hyperlogarithm.jl) and the Stokes constants by every method of `stokes_constant`
# (src/Mould.jl): the hyperlogarithm series truncated at weight N (`method = :hyperlog`, also with the
# certified choice N = :auto), the single-loop formula (`:loop`) and the matrix continuation (`:matrix`).
# The cost of the hyperlogarithm series grows as (letters)^N, so the tests stay at low weight.

const IR = IntegralResurgence

"Largest |coefficient| of the series A - B."
series_diff(A, B) = maximum(abs(get(A, w, 0) - get(B, w, 0)) for w in union(keys(A), keys(B)))

@testset "hyperlogarithms" begin
    # letters {0, 1}: L_0 = log z, L_1 = log(1 - z), L_{0 1} = -Li₂(z), shuffle L_0 L_1 = L_{01} + L_{10}
    Σ = ComplexF64[0, 1]
    z = 0.5 + 0.0im
    L = IR.base_value(Σ, z, 3)
    Li2(x) = sum(x^k / k^2 for k in 1:200)
    @test abs(L[[1]] - log(z)) < 1e-14
    @test abs(L[[2]] - log(1 - z)) < 1e-14
    @test abs(L[[1, 2]] + Li2(z)) < 1e-14
    @test abs(L[[1]] * L[[2]] - L[[1, 2]] - L[[2, 1]]) < 1e-14

    # one counterclockwise loop around z = 1: L_1 jumps by 2πi, L_0 does not
    Lγ = IR.continue_along(Σ, z, [(2, 1)], 3)
    @test abs(Lγ[[2]] - L[[2]] - 2π * im) < 1e-12
    @test abs(Lγ[[1]] - L[[1]]) < 1e-12
    # the monodromy factor C⁻¹ exp(2πi x₁) C agrees with explicit transport around the loop
    M, _ = IR.monodromy_factor(Σ, 2, z, 3)
    @test series_diff(M, IR.direct_loop_factor(Σ, 2, z, 3)) < 1e-12
    # a loop and its inverse cancel
    @test series_diff(IR.continue_along(Σ, z, [(2, 1), (2, -1)], 3), L) < 1e-12

    # series algebra
    A = IR.add!(IR.one_s(), IR.letter(1), 0.5); IR.add!(A, IR.letter(2), -0.25)
    @test series_diff(IR.mul(A, IR.sinv(A, 4), 4), IR.one_s()) < 1e-14
    x1 = IR.letter(1)
    @test series_diff(IR.mul(IR.sexp(x1, 4), IR.sexp(IR.scale(x1, -1), 4), 4), IR.one_s()) < 1e-14
    @test_throws ErrorException IR.transport(Σ, IR.one_s(), [0.5 + 0im, 1.0 + 0im], 2)   # runs into σ = 1
end

@testset "monodromy matrices from hyperlogarithms" begin
    # Y' = (R₁/z + R₂/(z - 1)) Y with nilpotent residues. For a nilpotent residue, the based monodromy
    # T_j = G exp(2πi R_j) G⁻¹ satisfies (T_j - I)² = 0 and T_j - I has rank one. The weight truncation
    # gains a factor of about 2π|R| per weight, so the residues are small.
    Σ = ComplexF64[0, 1]
    Rs = [0.01 * ComplexF64[0 1; 0 0], 0.01 * ComplexF64[0 0; 1 0]]
    z0 = 0.4 + 0.1im
    N = 4
    D = IR.monodromy_pieces(Σ, Rs, z0, N)
    T = IR.monodromy_matrices(Σ, Rs, z0, N)
    for j in 1:2
        @test T[j] ≈ sum(D[j])
        @test D[j][2] ≈ 2π * im * Rs[j]                       # weight one: 2πi R_j
        @test opnorm((T[j] - I)^2) < 1e-9 * opnorm(T[j] - I)^2
        s = svdvals(T[j] - I)
        @test s[2] < 1e-9 * s[1]                              # rank one
        # C⁻¹ exp(2πi x_j) C agrees with explicit transport around the loop, weight by weight
        @test series_diff(IR.monodromy_factor(Σ, j, z0, N)[1], IR.direct_loop_factor(Σ, j, z0, N)) < 1e-11
    end
    tc = IR.trace_contributions(D[1], D[2], N)
    @test abs(tc[1] - (2π * im)^2 * tr(Rs[1] * Rs[2])) < 1e-14
end

@testset "Stokes constants by every method (Pearcey integral)" begin setprecision(BigFloat, 256) do
    @polyvar q
    f = q^4 + 1//2 * q^2 + q
    s7 = sqrt(big(7))
    saddles = [[complex(-big(1)/2)], [(1 - s7 * im)/4], [(1 + s7 * im)/4]]
    GM = GaussManin(f)
    trans = transseries(f, GM, 10; saddles = saddles)
    Σ = ComplexF64[s.f₀ for s in trans]
    T0 = ComplexF64[s.Ts[1] for s in trans]
    Rs = [ComplexF64.(R) for R in GM.Rs]

    K = stokes_constants(Σ, Rs, 0; T0 = T0, method = :matrix)
    Kint = round.(Int, real.(K))
    @test maximum(abs, K - Kint) < 1e-10
    @test Kint == -transpose(Kint) && all(abs.(Kint[i, j]) == 1 for i in 1:3, j in 1:3 if i != j)

    # single-loop formula
    @test maximum(abs, stokes_constants(Σ, Rs, 0; T0 = T0, method = :loop) - K) < 1e-10

    # hyperlogarithm series at weight 4: within 1e-3 of the integers (the error decays with the weight)
    K4 = stokes_constants(Σ, Rs, 4; T0 = T0, method = :hyperlog)
    @test maximum(abs, K4 - K) < 1e-3
    K3 = stokes_constants(Σ, Rs, 3; T0 = T0, method = :hyperlog)
    @test maximum(abs, K4 - K) < maximum(abs, K3 - K)
    # the propagators of both back-ends agree along the keyhole path
    γp, _ = IR.keyhole_paths(Σ, 1, 2)
    Uh = IR.evaluate_series(IR.transport(Σ, IR.one_s(), γp, 4), Rs, 4)
    Um = IR.matrix_transport(Σ, Rs, γp)
    @test opnorm(Uh - Um) < 1e-2 * opnorm(Um)

    # certified truncation: B(N) ≥ |K_N - K| decreases with N; the first N with B(N) < 1/2 rounds exactly
    B = IR.truncation_bound(Σ, Rs, 1, 2; T0 = T0)
    @test B(2) > B(4) > B(6)
    @test abs(K4[2, 1] - K[2, 1]) <= B(4)
    Bprior = IR.truncation_bound(Σ, Rs, 1, 2; T0 = T0, a_posteriori = false, sharp = false)
    @test Bprior(6) > Bprior(10) > 0
    Nstar = IR.sufficient_N(Σ, Rs, 1, 2; T0 = T0)
    @test Nstar isa Int && B(Nstar) < 1/2
    Kauto = @test_logs (:info, r"certified") stokes_constant(Σ, Rs, 1, 2, :auto; T0 = T0)[1]
    @test round(Int, real(Kauto)) == Kint[2, 1]
    # word budget below the certified weight: falls back to the increment rule, with a warning
    Kbudget = @test_logs (:warn, r"exceeds the budget") match_mode = :any stokes_constant(Σ, Rs, 1, 2, :auto;
                                                                                    T0 = T0, max_words = 30)
    @test round(Int, real(Kbudget[1])) == Kint[2, 1]

    @test_throws ErrorException stokes_constant(Σ, Rs, 1, 2, 0; method = :nope)
end end
