using TensorKit, LinearAlgebra, MPSKit, KrylovKit
using JLD2, CairoMakie
using Revise
using IPEPSC6v

V = SU2Space(1//2=>1, 0=>1)
P = SU2Space(1//2=>1)
T = TensorMap(zeros, ComplexF64, P, V^4)

λs = 0.31:0.01:0.39
Etots = map(λs) do λ 
    Tfull, TA, TB, A, B = IPEPSC6v.long_range_RVB(λ)
    @load "data/long_range_RVB_lambda$(λ)_chi200.jld2" ψ1 ψ2
    # transfer matrix
    ψA = ψ2.AL[1]
    Etot, E1, E2, E3, E4 = IPEPSC6v.long_range_RVB_energy(Tfull, A, TB, ψA);
    return Etot
end
Etots_chi100 = map(λs) do λ 
    Tfull, TA, TB, A, B = IPEPSC6v.long_range_RVB(λ)
    @load "data/long_range_RVB_lambda$(λ)_chi100.jld2" ψ1 ψ2
    # transfer matrix
    ψA = ψ2.AL[1]
    Etot, E1, E2, E3, E4 = IPEPSC6v.long_range_RVB_energy(Tfull, A, TB, ψA);
    return Etot
end

fig = Figure(backgroundcolor=:white, fontsize=18, size=(600, 600))
ax1 = Axis(fig[1, 1], 
        xlabel = L"λ",
        ylabel = L"E", 
        )
lines!(ax1, λs, Etots, label="χ=200")
lines!(ax1, λs, Etots_chi100, label="χ=100")
axislegend(ax1, position=:rt)
@show fig

save("RVB_energies.pdf", fig)