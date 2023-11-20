module IPEPSC6v

__precompile__(true)

using LinearAlgebra
using TensorKit

export mapping_table, get_data_by_index, set_data_by_vector!, set_data_by_vector, selector, spatial_operation

# Write your package code here.
include("spatial_operations.jl");


end
