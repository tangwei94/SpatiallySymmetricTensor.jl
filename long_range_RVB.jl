# https://doi.org/10.1103/PhysRevLett.111.037202
# https://doi.org/10.1103/PhysRevB.94.205124

using TensorKit, LinearAlgebra, MPSKit, KrylovKit
using JLD2
using Revise
using IPEPSC6v

V = SU2Space(1//2=>1, 0=>1)
P = SU2Space(1//2=>1)
T = TensorMap(zeros, ComplexF64, P, V^4)

λ = 0.35#parse(Float64, ARGS[1])
for χ in 36:36:144
    Tfull, TA, TB = IPEPSC6v.long_range_RVB(λ)
    @show IPEPSC6v.mpo_hermicity(Tfull)
    @show IPEPSC6v.mpo_normality(Tfull)

    ψi = InfiniteMPS([fuse(V'*V)], [SU2Space(0=>1)])
    Tfull, TA, TB, A, B = IPEPSC6v.long_range_RVB(λ)
    𝕋full = DenseMPO([Tfull]) 
    let ψ1 = ψi, ψ2 = ψi
        ψ1 = ψi
        for ix in 1:100 
            ψ1 = changebonds(𝕋full * ψ1, SvdCut(truncdim(χ))) 
            @show ix, domain(ψ1.CR[1]) 
        end 
        ψ2, _, _ = leading_boundary(ψ1, 𝕋full, VUMPS(tol_galerkin=1e-12, maxiter=10000)); 
        @save "data/long_range_RVB_lambda$(λ)_chi$(χ).jld2" ψ1 ψ2
    end

    @load "data/long_range_RVB_lambda$(λ)_chi$(χ).jld2" ψ1 ψ2
    ## transfer matrix
    ψA = ψ2.AL[1]
    Etot, E1, E2, E3, E4 = IPEPSC6v.long_range_RVB_energy(Tfull, A, TB, ψA);
    @show Etot, E1, E2, E3, E4

    io = open("tmpdata.txt", "a");
    write(io, "$(λ) $(χ) $(E1) $(E2) $(E3) $(E4) $(Etot)\n")
    close(io)
end