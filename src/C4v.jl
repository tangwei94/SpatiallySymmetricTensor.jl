const C4v_σd_permutation = ((1, ), (3, 2, 5, 4))
const C4v_σv_permutation = ((1, ), (4, 3, 2, 5))
const C4v_R_permutation = ((1, ), (3, 4, 5, 2))

# in the order of σd, σv, R, see http://symmetry.jacobs-university.de/cgi-bin/group.cgi?group=404&option=4
const C4v_A1_reps = (1, 1, 1) 
const C4v_A2_reps = (-1, -1, 1) 
const C4v_B1_reps = (-1, 1, 1) 
const C4v_B2_reps = (1, -1, 1) 

function get_reps(name::Symbol)
    (name == :A1) && return C4v_A1_reps
    (name == :A2) && return C4v_A2_reps
    (name == :B1) && return C4v_B1_reps
    (name == :B2) && return C4v_B2_reps
end

function find_solution_C4v(T::AbstractTensorMap, reps_name::Symbol, P_filter=nothing)
    reps = get_reps(reps_name)

    mt = mapping_table(T)
    R_mat = spatial_operation(T, C4v_R_permutation; _mapping_table=mt)
    σd_mat = spatial_operation(T, C4v_σd_permutation; _mapping_table=mt)
    σv_mat = spatial_operation(T, C4v_σv_permutation; _mapping_table=mt)

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