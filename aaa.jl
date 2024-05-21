using TensorKit, LinearAlgebra, MPSKit, KrylovKit
using JLD2, CairoMakie
using Revise
using IPEPSC6v

V = SU2Space(1//2=>1, 0=>1);
P = SU2Space(1//2=>1);
T = TensorMap(zeros, ComplexF64, P, V^4);

for (ix, (f1, f2)) in enumerate(fusiontrees(T)) 
    if ix == 1
        #@show T[f1, f2], vec(T[f1, f2])
    end
        @show T[f1, f2], vec(T[f1, f2])
end
@show T

mapping_table(T::AbstractTensorMap)
new_mapping_table(T)

vec(T)

v = rand(12)
set_data_by_vector!(T, v)
vec(T) == v
        
condition(f1, f2) = length(findall(rep-> rep == SU2Irrep(1//2), f2.uncoupled)) == 1

selector(T, condition)