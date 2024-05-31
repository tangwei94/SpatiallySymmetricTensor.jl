using TensorKit, LinearAlgebra, MPSKit, KrylovKit
using JLD2, CairoMakie
using Revise
using SpatiallySymmetricTensor

V = SU2Space(1//2=>1, 0=>1)
P = SU2Space(1//2=>1)
T = TensorMap(zeros, ComplexF64, P, V^4)

λ = 0.35
χ = 144
Tfull, TA, TB, A, B = SpatiallySymmetricTensor.long_range_RVB(λ)
@load "data/long_range_RVB_lambda$(λ)_chi$(χ).jld2" ψ1 ψ2
# transfer matrix
ψA = ψ2.AL[1]
Etot, E1, E2, E3, E4 = SpatiallySymmetricTensor.long_range_RVB_energy(Tfull, A, TB, ψ2.AL[1]);
println("$(λ) $(Etot) $(E1) $(E2) $(E3) $(E4) ") 
Etot, E1, E2, E3, E4 = SpatiallySymmetricTensor.long_range_RVB_energy(Tfull, A, TB, ψ2.AR[1]);
println("$(λ) $(Etot) $(E1) $(E2) $(E3) $(E4) ") 


xs = [100, 200, 300]
ys = [-0.47822557886077166, -0.47825844598447487, -0.47826879344239914]

fig = Figure(backgroundcolor=:white, fontsize=18, size=(600, 600)); 
ax1 = Axis(fig[1, 1], xlabel = L"χ", ylabel=L"E")
lines!(ax1, 1 ./ xs, ys)
@show fig