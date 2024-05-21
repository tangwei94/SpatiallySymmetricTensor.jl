using TensorKit, LinearAlgebra, MPSKit, KrylovKit
using JLD2, CairoMakie
using Revise
using IPEPSC6v

V = SU2Space(1//2=>1, 0=>1);
P = SU2Space(1//2=>1);
T = TensorMap(zeros, ComplexF64, P, V^4);

fusiontrees(T)
for (ix, (f1, f2)) in enumerate(fusiontrees(T)) 
    if ix == 1
        T[f1, f2][1] = 1
    end
    @show ix, T[f1, f2]
end
@show T

mapping_table(T::AbstractTensorMap)
