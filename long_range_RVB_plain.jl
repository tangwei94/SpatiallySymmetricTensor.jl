# https://doi.org/10.1103/PhysRevLett.111.037202
# https://doi.org/10.1103/PhysRevB.94.205124

using TensorKit, LinearAlgebra, MPSKit, KrylovKit
using JLD2
using Revise
using IPEPSC6v

λ = 0.35#parse(Float64, ARGS[1])
Tfull_plain, TA_plain, TB_plain, A_plain, B_plain = IPEPSC6v.long_range_RVB(λ; use_symmetric_tensor=false)

let Tfull = Tfull_plain, A = A_plain, TB = TB_plain
    for χ in 36:36:144
        ψ2 = InfiniteMPS([ℂ^9], [ℂ^χ])
        𝕋full = DenseMPO([Tfull]) 
        if χ > 36
            @load "data/plain/long_range_RVB_lambda$(λ)_chi$(χ-36).jld2" ψ2
            ψ2, _ = changebonds(ψ2, 𝕋full, OptimalExpand(trscheme=truncdim(36)))
        end
        @show domain(ψ2.AL[1])
        ψ2, _, _ = leading_boundary(ψ2, 𝕋full, VUMPS(tol_galerkin=1e-12, maxiter=10000)); 
        @save "data/plain/long_range_RVB_lambda$(λ)_chi$(χ).jld2" ψ2

        @load "data/plain/long_range_RVB_lambda$(λ)_chi$(χ).jld2" ψ2
        ## transfer matrix
        ψA = ψ2.AL[1]
        Etot, E1, E2, E3, E4 = IPEPSC6v.long_range_RVB_energy(Tfull, A, TB, ψA; use_symmetric_tensor=false);
        @show Etot, E1, E2, E3, E4

        io = open("data/plaintmpdata.txt", "a");
        write(io, "$(λ) $(χ) $(E1) $(E2) $(E3) $(E4) $(Etot)\n")
        close(io)
    end
end

#Tfull, TA, TB, A, B = IPEPSC6v.long_range_RVB(λ; use_symmetric_tensor=true);
#let Tfull = Tfull, A = A, TB = TB  
#    @load "data/long_range_RVB_lambda$(λ)_chi$(36).jld2" ψ2
#    ψA = ψ2.AL[1]
#    Etot, E1, E2, E3, E4 = IPEPSC6v.long_range_RVB_energy(Tfull, A, TB, ψA; use_symmetric_tensor=true);
#    @show Etot, E1, E2, E3, E4
#end