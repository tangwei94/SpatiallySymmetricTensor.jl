# SpatiallySymmetricTensor

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://tangwei94.github.io/SpatiallySymmetricTensor.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://tangwei94.github.io/SpatiallySymmetricTensor.jl/dev/)
[![Build Status](https://github.com/tangwei94/SpatiallySymmetricTensor.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/tangwei94/SpatiallySymmetricTensor.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/tangwei94/SpatiallySymmetricTensor.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/tangwei94/SpatiallySymmetricTensor.jl)

This package aims to provide a simple way to combine point group symmetries (e.g., C4v, C3v) with internal symmetries (e.g., SU(2), U(1)) within the framework of [TensorKit.jl](https://github.com/Jutho/TensorKit.jl). 

A minimal example of using this package is as follows. 
Suppose we aim to generate a SU(2)-symmetric PEPS on the square lattice while ensuring that the local tensor remains invariant under reflections and rotations.
The internal SU(2) symmetry will be handled by TensorKit.jl, while the spatial symmetries can further be handled by this package:

```julia 
using TensorKit
using SpatiallySymmetricTensor 
```

We first specify the physical/virtual space for the PEPS tensor
```julia
# specify the virtual space of the PEPS tensor
V = SU2Space(1//2=>1, 0=>1) 
# specify the physical space of the PEPS tensor
P = SU2Space(1//2=>1) 
# generate a zero symmetric tensor
T0 = zeros(ComplexF64, P, V^4)  
```

Compared to a plain tensor, the SU(2)-symmetric tensor contains fewer free parameters.
Imposing spatial symmetries will further reduce the number of free parameters in the symmetric tensor.
In our case, we want to look for the symmetric tensors that belong to the A1 representation of the point group C4v. 
```julia
sols = find_solutions(C4v(), T0, :A1)
```
The returned `sols` is a vector of linear independent symmetric tensors which (i) have the same internal symmetric structure as `T0` and (ii) belong to the A1 representation of C4v. 
`sols` will contain only two entries, implying that we are left with only one free parameter after imposing the spatial symmetry. 
In the end, we can parametrize the PEPS tensor as, for example, 
```julia
function PEPS_tensor(x) 
    coefficients = [1.0, x]
    return sum(coefficients .* sol)
end
```
where `x` is the free parameter.