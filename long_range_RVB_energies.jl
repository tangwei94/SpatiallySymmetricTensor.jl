using TensorKit, LinearAlgebra, MPSKit, KrylovKit
using JLD2, CairoMakie
using Revise
using IPEPSC6v

V = SU2Space(1//2=>1, 0=>1)
P = SU2Space(1//2=>1)
T = TensorMap(zeros, ComplexF64, P, V^4)

λs = 0.02:0.02:0.26

function f_Etots(χ::Int)
    Etots = map(λs) do λ 
        Tfull, TA, TB, A, B = IPEPSC6v.long_range_RVB(λ)
        @load "data/experiment_with_symmetry_sectors_more_iTEBD_steps/long_range_RVB_lambda$(λ)_chi$(χ).jld2" ψ1 ψ2
        # transfer matrix
        ψA = ψ2.AL[1]
        Etot, E1, E2, E3, E4 = IPEPSC6v.long_range_RVB_energy(Tfull, A, TB, ψA);
        println("$(λ) $(Etot)   \r");
        return Etot
    end
    return Etots
end

Etots_chi100 = f_Etots(100)
Etots_chi200 = f_Etots(200)

fig = Figure(backgroundcolor=:white, fontsize=18, size=(600, 600))
ax1 = Axis(fig[1, 1], 
        xlabel = L"λ",
        ylabel = L"E", 
        )
lines!(ax1, λs, Etots_chi100, label="χ=100")
lines!(ax1, λs, Etots_chi200, label="χ=200")
axislegend(ax1, position=:rt)
@show fig

save("RVB_energies.pdf", fig)