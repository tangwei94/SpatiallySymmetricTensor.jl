abstract type AbstractSpaceGroup end

function find_solution(spg::Type{AbstractSpaceGroup}, T::AbstractTensorMap, reps_name::Symbol; P_filter=nothing)
    reps = get_reps(spg, reps_name)

    mt = mapping_table(T)
    R_mat = spatial_operation(T,  get_perm(spg, :R); _mapping_table=mt)
    σd_mat = spatial_operation(T, get_perm(spg, σd); _mapping_table=mt)
    σv_mat = spatial_operation(T, get_perm(spg, σv); _mapping_table=mt)

    if isnothing(P_filter)
        num_paras = num_free_parameters(T; _mapping_table=mt)
        P_filter = Matrix{ComplexF64}(I, num_paras, num_paras)
    end

    σd_mat1 = P_filter' * σd_mat * P_filter
    σv_mat1 = P_filter' * σv_mat * P_filter
    R_mat1 = P_filter' * R_mat * P_filter

    # σd
    Λσd, Uσd = eigen(Hermitian(σd_mat1))
    Pσd = Uσd[:, Λσd .≈ reps[1]]

    σv_mat2 = Pσd' * σv_mat1 * Pσd
    R_mat2 = Pσd' * R_mat1 * Pσd

    # σv
    Λσv, Uσv = eigen(Hermitian(σv_mat2))
    Pσv = Uσv[:, Λσv .≈ reps[2]]

    R_mat3 = Pσv' * R_mat2 * Pσv
    
    # R
    ΛR, UR = eigen(R_mat3)
    PR = UR[:, ΛR .≈ reps[3]]

    # solution 
    solutions = P_filter * Pσd * Pσv * PR
    num_solutions = size(solutions, 2)
    return [set_data_by_vector(T, vec(solutions[:, ix]); _mapping_table=mt) for ix in 1:num_solutions]
end