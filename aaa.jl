using TensorKit, LinearAlgebra, MPSKit, KrylovKit
using JLD2, CairoMakie
using Revise
using SpatiallySymmetricTensor
   
α = -1
T_3_1 = SpatiallySymmetricTensor.T_3_1_A1()
V = SU2Space(0=>1, 1//2=>1)
λ = TensorMap(zeros, ComplexF64, V, V); 
λ.data.values[1] .= 1
λ.data.values[2] .= α
@show convert(Array, λ)

@tensor A[-1; -2 -3 -4 -5] := T_3_1[-1; 1 2 3 4] * λ[1; -2] * λ[2; -3] * λ[3; -4] * λ[4; -5];
@tensor A2[-1; -2 -3 -4 -5] := T_3_1[-1; 1 2 -4 -5] * λ[1; -2] * λ[2; -3];
vec(A)
vec(T_3_1)

norm(A)/norm(T_3_1) * T_3_1 + A |> norm

norm(A)/norm(T_3_1) / α