using TensorKit, LinearAlgebra, MPSKit
using Revise
using IPEPSC6v

V = SU2Space(1//2=>1, 0=>1)
P = SU2Space(1//2=>1)
T = TensorMap(zeros, ComplexF64, P, V^4)

# a projector to the subspace with (1//2 ⊕ 0 ⊕ 0 ⊕ 0 -> 1//2)
P_nocc_1_3 = begin
    _condition(f1, f2) = length(findall(rep-> rep == SU2Irrep(1//2), f2.uncoupled)) == 1
    selector(T, _condition)
end
# a projector to the subspace with (1//2 ⊕ 1//2 ⊕ 1//2 ⊕ 0 -> 1//2)
P_nocc_3_1 = begin
    _condition(f1, f2) = length(findall(rep-> rep == SU2Irrep(1//2), f2.uncoupled)) == 3
    selector(T, _condition)
end

R_mat = spatial_operation(T, ((1, ), (3, 4, 5, 2)))
σd_mat = spatial_operation(T, ((1, ), (3, 2, 5, 4)))
σv_mat = spatial_operation(T, ((1, ), (4, 3, 2, 5)))

##### nocc= {1, 3}, TABLE VII in PRB 94, 205124 (2016)
T_1_3_A1 = begin
    Λσv_1_3, Uσv_1_3 = eigen(Hermitian(P_nocc_1_3' * σv_mat * P_nocc_1_3))
    Pσv_1_3 = Uσv_1_3[:, Λσv_1_3 .≈ 1]

    Λσd_1_3, Uσd_1_3 = eigen(Hermitian(Pσv_1_3' * P_nocc_1_3' * σd_mat * P_nocc_1_3 * Pσv_1_3))
    Pσd_1_3 = Uσd_1_3[:, Λσd_1_3 .≈ 1]

    ΛR_1_3, UR_1_3 = eigen(Pσd_1_3' * Pσv_1_3' * P_nocc_1_3' * R_mat * P_nocc_1_3 * Pσv_1_3 * Pσd_1_3)
    PR_1_3 = UR_1_3[:, ΛR_1_3 .≈ 1]

    sol_1_3 = vec(P_nocc_1_3 * Pσv_1_3 * Pσd_1_3 * PR_1_3)
    T_1_3_A1 = set_data_by_vector(T, sol_1_3)
end
##### nocc= {3, 1}, TABLE IX.  in PRB 94, 205124 (2016)
T_3_1_A1 = begin
    Λσv_3_1, Uσv_3_1 = eigen(Hermitian(P_nocc_3_1' * σv_mat * P_nocc_3_1))
    Pσv_3_1 = Uσv_3_1[:, Λσv_3_1 .≈ 1]

    Λσd_3_1, Uσd_3_1 = eigen(Hermitian(Pσv_3_1' * P_nocc_3_1' * σd_mat * P_nocc_3_1 * Pσv_3_1))
    Pσd_3_1 = Uσd_3_1[:, Λσd_3_1 .≈ 1]

    ΛR_3_1, UR_3_1 = eigen(Pσd_3_1' * Pσv_3_1' * P_nocc_3_1' * R_mat * P_nocc_3_1 * Pσv_3_1 * Pσd_3_1)
    PR_3_1 = UR_3_1[:, ΛR_3_1 .≈ 1]

    sol_3_1 = vec(P_nocc_3_1 * Pσv_3_1 * Pσd_3_1 * PR_3_1)
    T_3_1_A1 = set_data_by_vector(T, sol_3_1)
end

λ = 0.1
A = T_1_3_A1 + λ * T_3_1_A1
B = Tensor(zeros, ComplexF64, V*V)

@tensor TA[-1 -2 -3 -4; -5 -6 -7 -8] := A[1; -5 -6 -7 -8] * conj(A[1; -1 -2 -3 -4]) 
@tensor TB[-1 -2; -3 -4] := B[-1 -2] * conj(B[-3 -4])

δ = isomorphism(fuse(V'*V), V'*V);
@tensor Tfull[-1 -2; -3 -4] := TA[3 1 9 11; 4 2 10 12] * TB[6 2; 5 1] * TB[8 4; 7 3] * δ[-1; 5 6] * δ[-2; 7 8] * conj(δ[-3; 9 10]) * conj(δ[-4; 11 12]);
@show space(Tfull)

ψi = InfiniteMPS([fuse(V'*V)], [fuse(V'*V)])