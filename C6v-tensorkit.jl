using TensorKit, LinearAlgebra
using Revise
using IPEPSC6v

# example: (1//2 ⊕ 0)^6 -> 1//2
V = SU2Space(1//2=>1, 0=>1)
P = SU2Space(1//2=>1)
T = TensorMap(zeros, ComplexF64, P, V^6)

### projectors
P_nocc_1_5 = begin
    _condition(f1, f2) = length(findall(rep-> rep == SU2Irrep(1//2), f2.uncoupled)) == 1
    selector(T, _condition)
end
P_nocc_3_3 = begin
    _condition(f1, f2) = length(findall(rep-> rep == SU2Irrep(1//2), f2.uncoupled)) == 3
    selector(T, _condition)
end
P_nocc_5_1 = begin
    _condition(f1, f2) = length(findall(rep-> rep == SU2Irrep(1//2), f2.uncoupled)) == 5
    selector(T, _condition)
end

R_mat = spatial_operation(T, ((1, ), (3, 4, 5, 6, 7, 2)))
σd_mat = spatial_operation(T, ((1, ), (6, 5, 4, 3, 2, 7)))
σv_mat = spatial_operation(T, ((1, ), (7, 6, 5, 4, 3, 2)))

# IRREPS A1 : λR = λσd = λσv = 1
Ts_A1 = begin
    ##### nocc= {1, 5}
    Λσv_1_5, Uσv_1_5 = eigen(Hermitian(P_nocc_1_5' * σv_mat * P_nocc_1_5))
    Pσv_1_5 = Uσv_1_5[:, Λσv_1_5 .≈ 1]

    Λσd_1_5, Uσd_1_5 = eigen(Hermitian(Pσv_1_5' * P_nocc_1_5' * σd_mat * P_nocc_1_5 * Pσv_1_5))
    Pσd_1_5 = Uσd_1_5[:, Λσd_1_5 .≈ 1]

    ΛR_1_5, UR_1_5 = eigen(Pσd_1_5' * Pσv_1_5' * P_nocc_1_5' * R_mat * P_nocc_1_5 * Pσv_1_5 * Pσd_1_5)
    PR_1_5 = UR_1_5[:, ΛR_1_5 .≈ 1]

    sol_1_5 = vec(P_nocc_1_5 * Pσv_1_5 * Pσd_1_5 * PR_1_5)
    T_1_5 = set_data_by_vector(T, sol_1_5)
    for (f1, f2) in fusiontrees(T_1_5)
        (norm(T_1_5[f1, f2]) > 1e-12) && println(f1, " ", f2, " ", T_1_5[f1, f2])
    end

    ##### nocc= {3, 3}
    Λσv_3_3, Uσv_3_3 = eigen(Hermitian(P_nocc_3_3' * σv_mat * P_nocc_3_3))
    Pσv_3_3 = Uσv_3_3[:, Λσv_3_3 .≈ 1]

    Λσd_3_3, Uσd_3_3 = eigen(Hermitian(Pσv_3_3' * P_nocc_3_3' * σd_mat * P_nocc_3_3 * Pσv_3_3))
    Pσd_3_3 = Uσd_3_3[:, Λσd_3_3 .≈ 1]

    ΛR_3_3, UR_3_3 = eigen(Pσd_3_3' * Pσv_3_3' * P_nocc_3_3' * R_mat * P_nocc_3_3 * Pσv_3_3 * Pσd_3_3)
    PR_3_3 = UR_3_3[:, ΛR_3_3 .≈ 1]

    sol_3_3_a = vec((P_nocc_3_3 * Pσv_3_3 * Pσd_3_3 * PR_3_3)[:, 1])
    sol_3_3_b = vec((P_nocc_3_3 * Pσv_3_3 * Pσd_3_3 * PR_3_3)[:, 2])
    sol_3_3_c = vec((P_nocc_3_3 * Pσv_3_3 * Pσd_3_3 * PR_3_3)[:, 3])
    T_3_3_a = set_data_by_vector(T, sol_3_3_a)
    T_3_3_b = set_data_by_vector(T, sol_3_3_b)
    T_3_3_c = set_data_by_vector(T, sol_3_3_c)

    ##### nocc = {5, 1}
    Λσv_5_1, Uσv_5_1 = eigen(Hermitian(P_nocc_5_1' * σv_mat * P_nocc_5_1))
    Pσv_5_1 = Uσv_5_1[:, Λσv_5_1 .≈ 1]

    Λσd_5_1, Uσd_5_1 = eigen(Hermitian(Pσv_5_1' * P_nocc_5_1' * σd_mat * P_nocc_5_1 * Pσv_5_1))
    Pσd_5_1 = Uσd_5_1[:, Λσd_5_1 .≈ 1]

    ΛR_5_1, UR_5_1 = eigen(Pσd_5_1' * Pσv_5_1' * P_nocc_5_1' * R_mat * P_nocc_5_1 * Pσv_5_1 * Pσd_5_1)
    PR_5_1 = UR_5_1[:, ΛR_5_1 .≈ 1]

    sol_5_1_a = vec((P_nocc_5_1 * Pσv_5_1 * Pσd_5_1 * PR_5_1)[:, 1])
    sol_5_1_b = vec((P_nocc_5_1 * Pσv_5_1 * Pσd_5_1 * PR_5_1)[:, 2])
    sol_5_1_c = vec((P_nocc_5_1 * Pσv_5_1 * Pσd_5_1 * PR_5_1)[:, 3])
    T_5_1_a = set_data_by_vector(T, sol_5_1_a)
    T_5_1_b = set_data_by_vector(T, sol_5_1_b)
    T_5_1_c = set_data_by_vector(T, sol_5_1_c)
    Ts_A1 = [T_1_5, T_3_3_a, T_3_3_b, T_3_3_c, T_5_1_a, T_5_1_b, T_5_1_c]
end

# IRREPS A2 : λR = 1,  λσd = λσv = -1
Ts_A2 = begin
    Λσv_1_5, Uσv_1_5 = eigen(Hermitian(P_nocc_1_5' * σv_mat * P_nocc_1_5))
    Pσv_1_5 = Uσv_1_5[:, Λσv_1_5 .≈ -1]

    Λσd_1_5, Uσd_1_5 = eigen(Hermitian(Pσv_1_5' * P_nocc_1_5' * σd_mat * P_nocc_1_5 * Pσv_1_5))
    Pσd_1_5 = Uσd_1_5[:, Λσd_1_5 .≈ -1]
    # no solution

    ##### nocc= {3, 3}
    Λσv_3_3, Uσv_3_3 = eigen(Hermitian(P_nocc_3_3' * σv_mat * P_nocc_3_3))
    Pσv_3_3 = Uσv_3_3[:, Λσv_3_3 .≈ -1]

    Λσd_3_3, Uσd_3_3 = eigen(Hermitian(Pσv_3_3' * P_nocc_3_3' * σd_mat * P_nocc_3_3 * Pσv_3_3))
    Pσd_3_3 = Uσd_3_3[:, Λσd_3_3 .≈ -1]

    ΛR_3_3, UR_3_3 = eigen(Pσd_3_3' * Pσv_3_3' * P_nocc_3_3' * R_mat * P_nocc_3_3 * Pσv_3_3 * Pσd_3_3)
    PR_3_3 = UR_3_3[:, ΛR_3_3 .≈ 1]

    sol_3_3_a = vec((P_nocc_3_3 * Pσv_3_3 * Pσd_3_3 * PR_3_3)[:, 1])
    sol_3_3_b = vec((P_nocc_3_3 * Pσv_3_3 * Pσd_3_3 * PR_3_3)[:, 2])
    sol_3_3_c = vec((P_nocc_3_3 * Pσv_3_3 * Pσd_3_3 * PR_3_3)[:, 3])
    T_3_3_a = set_data_by_vector(T, sol_3_3_a)
    T_3_3_b = set_data_by_vector(T, sol_3_3_b)
    T_3_3_c = set_data_by_vector(T, sol_3_3_c)

    ##### nocc = {5, 1}
    Λσv_5_1, Uσv_5_1 = eigen(Hermitian(P_nocc_5_1' * σv_mat * P_nocc_5_1))
    Pσv_5_1 = Uσv_5_1[:, Λσv_5_1 .≈ -1]

    Λσd_5_1, Uσd_5_1 = eigen(Hermitian(Pσv_5_1' * P_nocc_5_1' * σd_mat * P_nocc_5_1 * Pσv_5_1))
    Pσd_5_1 = Uσd_5_1[:, Λσd_5_1 .≈ -1]

    ΛR_5_1, UR_5_1 = eigen(Pσd_5_1' * Pσv_5_1' * P_nocc_5_1' * R_mat * P_nocc_5_1 * Pσv_5_1 * Pσd_5_1)
    PR_5_1 = UR_5_1[:, ΛR_5_1 .≈ 1]

    sol_5_1_a = vec((P_nocc_5_1 * Pσv_5_1 * Pσd_5_1 * PR_5_1)[:, 1])
    sol_5_1_b = vec((P_nocc_5_1 * Pσv_5_1 * Pσd_5_1 * PR_5_1)[:, 2])
    T_5_1_a = set_data_by_vector(T, sol_5_1_a)
    T_5_1_b = set_data_by_vector(T, sol_5_1_b)
    
    Ts_A2 = [T_3_3_a, T_3_3_b, T_3_3_c, T_5_1_a, T_5_1_b]
end

λ1s = fill(1.0, 7)
T_A1 = sum(λ1s .* Ts_A1)
λ2s = fill(1.0, 5)
T_A2 = sum(λ2s .* Ts_A2)
T = T_A1 + im * T_A2