using TensorKit, LinearAlgebra, MPSKit, KrylovKit
using JLD2, CairoMakie
using Revise
using IPEPSC6v

V = SU2Space(1//2=>1, 0=>1)
P = SU2Space(1//2=>1)
T = TensorMap(zeros, ComplexF64, P, V^4)

λ = 0.34
χ = 100
Tfull, TA, TB, A, B = IPEPSC6v.long_range_RVB(λ)
@load "data/experiment_with_symmetry_sectors_more_iTEBD_steps/long_range_RVB_lambda$(λ)_chi$(χ).jld2" ψ1 ψ2
# transfer matrix
ψA = ψ2.AL[1]
Etot, E1, E2, E3, E4 = IPEPSC6v.long_range_RVB_energy(Tfull, A, TB, ψ1.AL[1]);
println("$(λ) $(Etot) $(E1) $(E2) $(E3) $(E4) ") 
Etot, E1, E2, E3, E4 = IPEPSC6v.long_range_RVB_energy(Tfull, A, TB, ψ2.AR[1]);
println("$(λ) $(Etot) $(E1) $(E2) $(E3) $(E4) ") 
