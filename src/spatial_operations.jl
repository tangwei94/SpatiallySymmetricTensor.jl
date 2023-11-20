function mapping_table(T::AbstractTensorMap)
    res = Tuple{Int, Int}[]
    for (ix, v) in enumerate(T.data.values)
        for iy in 1:length(v)
            push!(res, (ix, iy))
        end
    end
    return res
end

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

function selector(T::AbstractTensorMap, _condition; _mapping_table=mapping_table(T))
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