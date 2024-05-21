const MappingEntry = Tuple{FusionTree, FusionTree, Int, Int}
const MappingTable = Vector{MappingEntry}

"""
    mapping_table(T::AbstractTensorMap)

    - generate a mapping table for the free parameters in a symmetric tensor. Constructs the mapping between the linear space of the free parameters and the symmetric tensor.
"""
function mapping_table(T::AbstractTensorMap)
    res = MappingEntry[]
    a = 1
    for (f1, f2) in fusiontrees(T)
        n = length(T[f1, f2])
        push!(res, (f1, f2, a, n))
        a += n
    end
    return res
end

"""
    num_free_parameters(T::AbstractTensorMap)

    - return the number of free parameters in a symmetric tensor
"""
function num_free_parameters(T::AbstractTensorMap; _mapping_table::MappingTable=mapping_table(T))
    return sum([n for (_, _, _, n) in _mapping_table(T)])
end

"""
    set_data_by_vector!(T::AbstractTensorMap, _values; _mapping_table=mapping_table(T))

    - modify the free parameters in a symmetric tensor by values
    - map a vector that lives the linear space of the free parameters to a symmetric tensor. The mapping is defined by the symmetric structure of the symmetric tensor 
"""
function set_data_by_vector!(T::AbstractTensorMap, values::Vector{<:Number}; _mapping_table::MappingTable=mapping_table(T))
    for ix in eachindex(_mapping_table)
        f1, f2, a, n = _mapping_table[ix]
        T[f1, f2] .= values[a:a+n-1]
    end
    return T
end

"""
    set_data_by_vector(T::AbstractTensorMap, values::Vector{<:Number}; _mapping_table=mapping_table(T))

    - generate a new tensor which has the same symmetry structure as T, but with free parameters set by `values`
    - map a vector that lives the linear space of the free parameters to a symmetric tensor. The mapping is defined by the symmetric structure of the symmetric tensor 
"""
function set_data_by_vector(T::AbstractTensorMap, values::Vector{<:Number}; _mapping_table::MappingTable=mapping_table(T))
    T1 = similar(T) 
    set_data_by_vector!(T1, values; _mapping_table=_mapping_table)
    return T1
end

"""
    Base.vec(T::AbstractTensorMap; _mapping_table=mapping_table(T))

    - convert the free parameters in a tensor to a vector, in the same order as set_data_by_vector

    - map the symmetric tensor to a vector that lives the linear space of the free parameters
"""
function Base.vec(T::AbstractTensorMap; _mapping_table::MappingTable=mapping_table(T))
    v = map(eachindex(_mapping_table)) do ix
        f1, f2, _, _ = _mapping_table[ix]
        return vec(T[f1, f2])
    end
    return vcat(v...)
end

"""
    selector(T::AbstractTensorMap, _condition; _mapping_table=mapping_table(T))

    - construct a projector in the linear space of the free parameters of T
    - the subspace corresponds to the symmetric tensors that satisfy condition
returreturn
    - `_condition` is a function that takes two arguments `f1` and `f2`, which are the `FusionTrees` of `T`. The function return boolean values. For the explanation of `TensorKit.FusionTrees`, see https://jutho.github.io/TensorKit.jl/latest/man/sectors/.  
"""
function selector(T::AbstractTensorMap, condition::Function; _mapping_table::MappingTable=mapping_table(T))
    num_paras = sum([n for (_, _, _, n) in _mapping_table])
    P = zeros(ComplexF64, num_paras, num_paras)
    Pd = view(P, diagind(P))

    for ix in eachindex(_mapping_table)
        f1, f2, a, n = _mapping_table[ix]
        condition(f1, f2) && (Pd[a:a+n-1] .= 1)
    end
    return P
end

"""
    spatial_operation(T::AbstractTensorMap, permutations; _mapping_table=mapping_table(T))

    - When you do a spatial operation to a symmetric tensor, you map a symmetric tensor to another symmetric tensor. 
    - This mapping can be represented by a matrix defined in the linear space of the free parameters of the symmetric tensor.
    - The function `spatial_operation` returns this matrix.

"""
function spatial_operation(T::AbstractTensorMap, permutations; _mapping_table::MappingTable=mapping_table(T))
    num_paras = sum([n for (_, _, _, n) in _mapping_table])
    M = zeros(ComplexF64, num_paras, num_paras)
    for ix in 1:num_paras
        paras = zeros(ComplexF64, num_paras)
        paras[ix] = 1
    
        T1 = set_data_by_vector(T, paras; _mapping_table=_mapping_table)

        RT = permute(T1, permutations)
        X = vec(RT)
        M[:, ix] = X 
    end
    return M
end