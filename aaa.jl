using TensorKit, LinearAlgebra, MPSKit, KrylovKit
using JLD2, CairoMakie
using Revise
using IPEPSC6v
    
T13 = IPEPSC6v.T_1_3_A1() 
T13
T13_arr = convert(Array, T13)

λ13 = 0.5 / T13_arr[1, 2, 1, 1, 1] 
λ13*T13_arr[1, 2, 1, 1, 1]
λ13*T13_arr[1, 1, 2, 1, 1]
λ13*T13_arr[1, 1, 1, 2, 1]
λ13*T13_arr[1, 1, 1, 1, 2]
λ13*T13_arr[2, 3, 1, 1, 1]
λ13*T13_arr[2, 1, 3, 1, 1]
λ13*T13_arr[2, 1, 1, 3, 1]
λ13*T13_arr[2, 1, 1, 1, 3]

T31 = IPEPSC6v.T_3_1_A1()
T31_arr = convert(Array, T31)

λ31 = -(1/2/sqrt(6)) / T31_arr[1, 2, 2, 1, 3]
@show 1/sqrt(6)
@show 1/2/sqrt(6)

λ31*T31_arr[1, 2, 2, 1, 3]
λ31*T31_arr[1, 2, 2, 3, 1]
λ31*T31_arr[1, 2, 1, 2, 3]
λ31*T31_arr[1, 3, 2, 1, 2]

λ31*T31_arr[1, 3, 2, 1, 2]
λ31*T31_arr[1, 3, 2, 1, 2]
λ31*T31_arr[1, 3, 2, 1, 2]
λ31*T31_arr[1, 3, 2, 1, 2]

c = 0.35
λ13 * T13 + c * λ31 * T31

c * λ31 / λ13
