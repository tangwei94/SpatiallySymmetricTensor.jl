abstract type AbstractPointGroup end

function find_solution(spg::AbstractPointGroup, T::AbstractTensorMap, reps_name::Symbol; P_filter=nothing)
    reps = get_reps(spg, reps_name)

    mt = mapping_table(T)
    σd_mats = [spatial_operation(T, perm; _mapping_table=mt) for perm in get_perm(spg, :σd)] 
    σv_mats = [spatial_operation(T, perm; _mapping_table=mt) for perm in get_perm(spg, :σv)] 
    R_mats = [spatial_operation(T, perm; _mapping_table=mt) for perm in get_perm(spg, :R)] 

    if isnothing(P_filter)
        num_paras = num_free_parameters(T; _mapping_table=mt)
        P_filter = Matrix{ComplexF64}(I, num_paras, num_paras)
    end

    P_sol = P_filter

    for σd_mat in σd_mats
        (size(P_sol, 2) == 0) && break
        σd_mat1 = P_sol' * σd_mat * P_sol
        Λσd, Uσd = eigen(Hermitian(σd_mat1))
        Pσd = Uσd[:, Λσd .≈ reps[1]]
        P_sol = P_sol * Pσd
    end
    for σv_mat in σv_mats
        (size(P_sol, 2) == 0) && break
        σv_mat1 = P_sol' * σv_mat * P_sol
        Λσv, Uσv = eigen(Hermitian(σv_mat1))
        Pσv = Uσv[:, Λσv .≈ reps[2]]
        P_sol = P_sol * Pσv
    end
    for R_mat in R_mats
        (size(P_sol, 2) == 0) && break
        R_mat1 = P_sol' * R_mat * P_sol # R is not Hermitian
        ΛR, UR = eigen(R_mat1)
        PR = UR[:, ΛR .≈ reps[3]]
        P_sol = P_sol * PR
    end

    # solution 
    num_solutions = size(P_sol, 2)
    return [set_data_by_vector(T, vec(P_sol[:, ix]); _mapping_table=mt) for ix in 1:num_solutions]
end