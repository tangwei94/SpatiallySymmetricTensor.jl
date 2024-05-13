"""
    mapping_table(T::AbstractTensorMap)

    - generate a mapping table for the free parameters in a symmetric tensor. Constructs the mapping between the linear space of the free parameters and the symmetric tensor.
"""
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

"""
    set_data_by_vector!(T::AbstractTensorMap, _values; _mapping_table=mapping_table(T))

    - modify the free parameters in a symmetric tensor by _values
    - map a vector that lives the linear space of the free parameters to a symmetric tensor. The mapping is defined by the symmetric structure of the symmetric tensor 
"""
function set_data_by_vector!(T::AbstractTensorMap, _values; _mapping_table=mapping_table(T))
    for ix in eachindex(_values)
        T.data.values[_mapping_table[ix][1]][_mapping_table[ix][2]] = _values[ix]
    end
    T
end

"""
    set_data_by_vector(T::AbstractTensorMap, _values; _mapping_table=mapping_table(T))

    - generate a new tensor which has the same symmetry structure as T, but with free parameters set by `_values`
    - map a vector that lives the linear space of the free parameters to a symmetric tensor. The mapping is defined by the symmetric structure of the symmetric tensor 
"""
function set_data_by_vector(T::AbstractTensorMap, _values; _mapping_table=mapping_table(T))
    T1 = zero(T) 
    set_data_by_vector!(T1, _values; _mapping_table=_mapping_table)
end

"""
    Base.vec(T::AbstractTensorMap; _mapping_table=mapping_table(T))

    - convert the free parameters in a tensor to a vector, in the same order as set_data_by_vector

    - map the symmetric tensor to a vector that lives the linear space of the free parameters
"""
function Base.vec(T::AbstractTensorMap; _mapping_table=mapping_table(T))
    map(eachindex(_mapping_table)) do ix
        T.data.values[_mapping_table[ix][1]][_mapping_table[ix][2]]
    end
end

"""
    selector(T::AbstractTensorMap, _condition; _mapping_table=mapping_table(T))

    - construct a projector in the linear space of the free parameters of T
    - the subspace corresponds to the symmetric tensors that satisfy _condition
returreturnn
    - `_condition` is a function that takes two arguments `f1` and `f2`, which are the `FusionTrees` of `T`. The function return boolean values. For the explanation of `TensorKit.FusionTrees`, see https://jutho.github.io/TensorKit.jl/latest/man/sectors/.  
"""
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

"""
    spatial_operation(T::AbstractTensorMap, permutations; _mapping_table=mapping_table(T))

    - When you do a spatial operation to a symmetric tensor, you map a symmetric tensor to another symmetric tensor. 
    - This mapping can be represented by a matrix defined in the linear space of the free parameters of the symmetric tensor.
    - The function `spatial_operation` returns this matrix.

"""
function spatial_operation(T::AbstractTensorMap, permutations; _mapping_table=mapping_table(T))
    num_paras = length(_mapping_table)
    M = zeros(ComplexF64, num_paras, num_paras)
    for ix in 1:num_paras
        paras = zeros(ComplexF64, num_paras)
        paras[ix] = 1
    
        T1 = set_data_by_vector(T, paras; _mapping_table=_mapping_table)

        RT = permute(T1, permutations...)
        @assert fusiontrees(RT) == fusiontrees(T1)
        X = vec(RT)
        M[:, ix] = X 
    end
    return M
end