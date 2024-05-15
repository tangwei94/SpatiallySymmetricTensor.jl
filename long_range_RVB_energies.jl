using TensorKit, LinearAlgebra, MPSKit, KrylovKit
using JLD2, CairoMakie
using Revise
using IPEPSC6v

V = SU2Space(1//2=>1, 0=>1)
P = SU2Space(1//2=>1)
T = TensorMap(zeros, ComplexF64, P, V^4)

λ = 0.14
χs = 120:20:200
Etots = map(χs) do χ 
    Tfull, TA, TB, A, B = IPEPSC6v.long_range_RVB(λ)
    @load "data/itebd1_long_range_RVB_lambda$(λ)_chi$(χ).jld2" ψ1 ψ2
    # transfer matrix
    ψA = ψ2.AL[1]
    Etot, E1, E2, E3, E4 = IPEPSC6v.long_range_RVB_energy(Tfull, A, TB, ψA);
    return Etot
end

fig = Figure(backgroundcolor=:white, fontsize=18, size=(600, 600))
ax1 = Axis(fig[1, 1], 
        xlabel = L"χ",
        ylabel = L"E", 
        )
lines!(ax1, χs, Etots, label="χ=200")
axislegend(ax1, position=:rt)
@show fig

save("RVB_energies.pdf", fig)