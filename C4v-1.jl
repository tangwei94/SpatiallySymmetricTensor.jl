using TensorKit

A = zeros(3,2,2,2,2)
A[1,1,1,1,2]=1/2
A[1,1,1,2,1]=-1/2
A[1,1,2,1,1]=1/2
A[1,2,1,1,1]=-1/2

A[2,1,2,1,2]=1/√2
A[2,2,1,2,1]=-1/√2

A[3,1,2,2,2]=1/2
A[3,2,1,2,2]=-1/2
A[3,2,2,1,2]=1/2
A[3,2,2,2,1]=-1/2

V = SU2Space(1//2=>1)
P = SU2Space(1=>1)
PEPS = TensorMap(A, P, V^4) # SU2 symmetric
convert(Array, PEPS) # Convert back to normal array

#########################################

V = SU2Space(1//2=>1)
P = SU2Space(1//2=>1)
PEPS = TensorMap(zeros, ComplexF64, P, V^3) # SU2 symmetric
PEPS.data.keys
PEPS.data.values[1][1] = 1
PEPS.data

function mapping_table(T::AbstractTensorMap)
    res = Tuple{Int, Int}[]
    for (ix, v) in enumerate(T.data.values)
        for iy in 1:length(v)
            push!(res, (ix, iy))
        end
    end
    return res
end
tab = mapping_table(PEPS)

function get_data_by_index(T::AbstractTensorMap, ix::Int; _mapping_table=mapping_table(T))
    T.data.values[_mapping_table[ix][1]][_mapping_table[ix][2]]
end

# modify the free parameters in a tensor
function set_data_by_vector!(T::AbstractTensorMap, _values; _mapping_table=mapping_table(T))
    for ix in eachindex(_values)
        T.data.values[_mapping_table[ix][1]][_mapping_table[ix][2]] = _values[ix]
    end
    T
end
# generate a new tensor which has the same symmetry structure as T, but with free parameters set by _values
function set_data_by_vector(T::AbstractTensorMap, _values; _mapping_table=mapping_table(T))
    T1 = zero(T) 
    set_data_by_vector!(T1, _values; _mapping_table=_mapping_table)
end
# convert the free parameters in a tensor to a vector, in the same order as set_data_by_vector
function Base.vec(T::AbstractTensorMap; _mapping_table=mapping_table(T))
    map(eachindex(_mapping_table)) do ix
        T.data.values[_mapping_table[ix][1]][_mapping_table[ix][2]]
    end
end

## example: (1//2 ⊕ 0)^4 -> 1//2
V = SU2Space(1//2=>1, 0=>1) 
P = SU2Space(1//2=>1)
T = TensorMap(zeros, ComplexF64, P, V^4)
T

function filter(T::AbstractTensorMap, _condition; _mapping_table=mapping_table(T))
    num_paras = length(_mapping_table)
    P = zeros(ComplexF64, num_paras, num_paras)
    for ix in 1:num_paras
        paras = zeros(ComplexF64, num_paras)
        paras[ix] = 1
        T1 = set_data_by_vector(T, paras; _mapping_table=_mapping_table)
        
        for (f1, f2) in fusiontrees(T1)
            (norm(T1[f1, f2]) > 1e-12) && _condition(f1, f2) && (P[ix, ix] = 1)
        end
    end
    return P[:, findall(ix->P[ix, ix]==1, 1:num_paras)]
end

function spatial_operation(T::AbstractTensorMap, permutations; _mapping_table=mapping_table(T))
    num_paras = length(_mapping_table)
    M = zeros(ComplexF64, num_paras, num_paras)
    for ix in 1:num_paras
        paras = zeros(ComplexF64, num_paras)
        paras[ix] = 1
    
        T1 = set_data_by_vector(T, paras; _mapping_table=_mapping_table)
        @show T1

        RT = permute(T1, permutations...)
        @assert fusiontrees(RT) == fusiontrees(T1)
        X = vec(RT)
        M[:, ix] = X 
    end
    return M
end

# a projector to the subspace with (1//2 ⊕ 0 ⊕ 0 ⊕ 0 -> 1//2)
P_nocc_1_3 = begin
    _condition(f1, f2) = length(findall(rep-> rep == SU2Irrep(1//2), f2.uncoupled)) == 1
    filter(T, _condition)
end
# a projector to the subspace with (1//2 ⊕ 1//2 ⊕ 1//2 ⊕ 0 -> 1//2)
P_nocc_3_1 = begin
    _condition(f1, f2) = length(findall(rep-> rep == SU2Irrep(1//2), f2.uncoupled)) == 3
    filter(T, _condition)
end

R_mat = spatial_operation(T, ((1, ), (3, 4, 5, 2)))
σd_mat = spatial_operation(T, ((1, ), (3, 2, 5, 4)))
σv_mat = spatial_operation(T, ((1, ), (4, 3, 2, 5)))

# test
using Test
num_paras = length(mapping_table(T))
for ix0 in 1:num_paras
    T0 = begin
        paras = zeros(num_paras)
        paras[ix0] = 1
        set_data_by_vector(T, paras)
    end

    T1 = set_data_by_vector(T, R_mat[:, ix0])
    @test norm(permute(T0, (1, ), (3, 4, 5, 2)) - T1) < 1e-12
end


# IRREPS A1 : λR = λσd = 1 = λσv = 1

### does not use the projectors P_nocc_1_3 and P_nocc_3_1, got two variational parameters, but the result is ``mixed''
Λσv, Uσv = eigen(Hermitian(σv_mat))
Pσv = Uσv[:, Λσv .≈ 1]

Λσd, Uσd = eigen(Hermitian(Pσv' * σd_mat * Pσv))
Pσd = Uσd[:, Λσd .≈ 1]

ΛR, UR = eigen(Pσd' * Pσv' * R_mat * Pσv * Pσd)
PR = UR[:, ΛR .≈ 1]

sol1 = vec(Pσv * Pσd * PR[:, 1])
sol2 = vec(Pσv * Pσd * PR[:, 2])

PEPS1 = set_data_by_vector(T, sol1)
for (f1, f2) in fusiontrees(PEPS1)
    (norm(PEPS1[f1, f2]) > 1e-12) && println(f1, " ", f2, " ", PEPS1[f1, f2])
end

PEPS2 = set_data_by_vector(T, sol2)
for (f1, f2) in fusiontrees(PEPS2)
    (norm(PEPS2[f1, f2]) > 1e-12) && println(f1, " ", f2, " ", PEPS2[f1, f2])
end

### use the projectors P_nocc_1_3 and P_nocc_3_1
##### nocc= {1, 3}, TABLE VII in PRB 94, 205124 (2016)
Λσv_1_3, Uσv_1_3 = eigen(Hermitian(P_nocc_1_3' * σv_mat * P_nocc_1_3))
Pσv_1_3 = Uσv_1_3[:, Λσv_1_3 .≈ 1]

Λσd_1_3, Uσd_1_3 = eigen(Hermitian(Pσv_1_3' * P_nocc_1_3' * σd_mat * P_nocc_1_3 * Pσv_1_3))
Pσd_1_3 = Uσd_1_3[:, Λσd_1_3 .≈ 1]

ΛR_1_3, UR_1_3 = eigen(Pσd_1_3' * Pσv_1_3' * P_nocc_1_3' * R_mat * P_nocc_1_3 * Pσv_1_3 * Pσd_1_3)
PR_1_3 = UR_1_3[:, ΛR_1_3 .≈ 1]

sol_1_3 = vec(P_nocc_1_3 * Pσv_1_3 * Pσd_1_3 * PR_1_3)
PEPS_1_3 = set_data_by_vector(T, sol_1_3)
for (f1, f2) in fusiontrees(PEPS_1_3)
    (norm(PEPS_1_3[f1, f2]) > 1e-12) && println(f1, " ", f2, " ", PEPS_1_3[f1, f2])
end

##### nocc= {3, 1}, TABLE IX.  in PRB 94, 205124 (2016)
Λσv_3_1, Uσv_3_1 = eigen(Hermitian(P_nocc_3_1' * σv_mat * P_nocc_3_1))
Pσv_3_1 = Uσv_3_1[:, Λσv_3_1 .≈ 1]

Λσd_3_1, Uσd_3_1 = eigen(Hermitian(Pσv_3_1' * P_nocc_3_1' * σd_mat * P_nocc_3_1 * Pσv_3_1))
Pσd_3_1 = Uσd_3_1[:, Λσd_3_1 .≈ 1]

ΛR_3_1, UR_3_1 = eigen(Pσd_3_1' * Pσv_3_1' * P_nocc_3_1' * R_mat * P_nocc_3_1 * Pσv_3_1 * Pσd_3_1)
PR_3_1 = UR_3_1[:, ΛR_3_1 .≈ 1]

sol_3_1 = vec(P_nocc_3_1 * Pσv_3_1 * Pσd_3_1 * PR_3_1)
PEPS_3_1 = set_data_by_vector(T, sol_3_1)
for (f1, f2) in fusiontrees(PEPS_3_1)
    (norm(PEPS_3_1[f1, f2]) > 1e-12) && println(f1, " ", f2, " ", PEPS_3_1[f1, f2])
end