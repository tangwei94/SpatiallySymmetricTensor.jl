@testset "spatial_operations.jl" begin
    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)
    T = TensorMap(zeros, ComplexF64, P, V^6)
    R_mat = spatial_operation(T, ((1, ), (3, 4, 5, 6, 7, 2)))
    num_paras = length(mapping_table(T))
    for ix0 in 1:num_paras
        T0 = begin
            paras = zeros(num_paras)
            paras[ix0] = 1
            set_data_by_vector(T, paras)
        end

        T1 = set_data_by_vector(T, R_mat[:, ix0])
        @test norm(permute(T0, (1, ), (3, 4, 5, 6, 7, 2)) - T1) < 1e-12
    end
end