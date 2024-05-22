@testset "test set_data_by_vector and vec" for _ in 1:10 
    # TODO. more test cases
    V = SU2Space(1//2=>1, 0=>1);
    P = SU2Space(1//2=>1);
    T = TensorMap(zeros, ComplexF64, P, V^4);

    mt = mapping_table(T)
    num_paras = num_free_parameters(T; _mapping_table=mt)

    v = rand(num_paras)

    T1 = set_data_by_vector(T, v; _mapping_table=mt)
    @test norm(vec(T1) - v) < 1e-12
end

@testset "spatial_operations.jl" for _ in 1:10
    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)
    T = TensorMap(zeros, ComplexF64, P, V^6)

    permutations = [((1, ), (3, 4, 5, 6, 7, 2)), 
                    ((1, ), (4, 5, 6, 7, 2, 3)),
                    ((1, ), (5, 6, 7, 2, 3, 4)),
                    ((1, ), (6, 7, 2, 3, 4, 5)),
                    ((1, ), (7, 2, 3, 4, 5, 6)),
                    ((1, ), (4, 3, 6, 2, 7, 5))]

    for perm in permutations
        R_mat = spatial_operation(T, perm)
        num_paras = num_free_parameters(T)
        for ix0 in 1:num_paras
            T0 = begin
                paras = zeros(num_paras)
                paras[ix0] = 1
                set_data_by_vector(T, paras)
            end
        
            T1 = set_data_by_vector(T, R_mat[:, ix0])
            @test norm(permute(T0, perm) - T1) < 1e-12
        end
    end
end