using TensorKit, LinearAlgebra, MPSKit, KrylovKit
using JLD2, CairoMakie
using Revise
using IPEPSC6v
    
V = SU2Space(1//2=>1, 0=>1)
P = SU2Space(1//2=>1)
T = TensorMap(zeros, ComplexF64, P, V^4)
sols = find_solution(C4v(), T, :A1);
length(sols)

V = U1Space(-1//2=>1, -0=>3, 1//2=>1)
P = U1Space(0=>1, 1=>1)
T = TensorMap(zeros, ComplexF64, P, V^4)
num_free_parameters(T)
sols = find_solution(C4v(), T, :A1);
length(sols)

a = convert(Array, sols[1])
a - permutedims(a, (1, 4, 3, 2, 5)) |> norm
a - permutedims(a, (1, 2, 4, 3, 5)) |> norm
a - permutedims(a, (1, 5, 3, 4, 2)) |> norm
a - permutedims(a, (1, 3, 2, 5, 4)) |> norm
a - permutedims(a, (1, 3, 4, 5, 2)) |> norm



# a projector to the subspace with (1//2 ⊕ 0 ⊕ 0 ⊕ 0 -> 1//2) (short-range RVB)
P_nocc_1_3 = begin
    _condition(f1, f2) = length(findall(rep-> rep == SU2Irrep(1//2), f2.uncoupled)) == 1
    selector(T, _condition)
end

P_trivial = Matrix{ComplexF64}(I, num_paras, num_paras)

find_sol(T, P_trivial, C4v_A1_reps)

V = SU2Space(1//2=>1, 0=>1)
P = SU2Space(1//2=>1)
T = TensorMap(zeros, ComplexF64, P, V^4)



