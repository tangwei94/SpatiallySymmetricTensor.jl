# https://doi.org/10.1103/PhysRevLett.111.037202
# https://doi.org/10.1103/PhysRevB.94.205124

using TensorKit, LinearAlgebra, MPSKit, KrylovKit
using JLD2
using Revise
using IPEPSC6v

V = SU2Space(1//2=>1, 0=>1)
P = SU2Space(1//2=>1)
T = TensorMap(zeros, ComplexF64, P, V^4)

T_1_3_A1 = IPEPSC6v.T_1_3_A1()
T_3_1_A1 = IPEPSC6v.T_1_3_A1()

λ = parse(Float64, ARGS[1])
λ = 0.35
Tfull, TA, TB = IPEPSC6v.long_range_RVB(λ)
@show IPEPSC6v.mpo_hermicity(Tfull)
@show IPEPSC6v.mpo_normality(Tfull)

ψi = InfiniteMPS([fuse(V'*V)], [fuse(V'*V)])

Tfull, TA, TB, A, B = IPEPSC6v.long_range_RVB(λ)

𝕋full = DenseMPO([Tfull]) 
let ψ1 = ψi, ψ2 = ψi
    ψ1 = ψi
    for ix in 1:20 
        ψ1 = changebonds(𝕋full * ψ1, SvdCut(truncdim(200))) 
        @show ix, domain(ψ1.CR[1]) 
    end 
    ψ2, _, _ = leading_boundary(ψ1, 𝕋full, VUMPS(tol_galerkin=1e-12, maxiter=1000)); 
    @save "data/long_range_RVB_lambda$(λ)_chi200.jld2" ψ1 ψ2
end

@load "data/long_range_RVB_lambda$(λ)_chi200.jld2" ψ1 ψ2
# transfer matrix
ψA = ψ2.AL[1]
Etot, E1, E2, E3, E4 = IPEPSC6v.long_range_RVB_energy(Tfull, A, TB, ψA);
@show Etot, E1, E2, E3, E4

io = open("tmpdata.txt", "a");
write(io, "$(λ) $(E1) $(E2) $(E3) $(E4) $(Etot)\n")
close(io)