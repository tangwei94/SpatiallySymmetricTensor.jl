function mpo_ovlp(A1, A2)
    V1 = MPSKit._firstspace(A1)
    V2 = MPSKit._firstspace(A2)

    function mpo_transf(v)
        @tensor Tv[-1; -2] := A1[-1 3; 4 1] * conj(A2[-2 3; 4 2]) * v[1; 2]
        return Tv
    end

    @show V1, V2
    v0 = TensorMap(rand, ComplexF64, V1, V2)
    return eigsolve(mpo_transf, v0, 1, :LM)
end

function mpotensor_dag(T::MPSKit.MPOTensor)
    T_data = reshape(T.data, (dims(codomain(T))..., dims(domain(T))...))
    Tdag_data = permutedims(conj.(T_data), (1, 3, 2, 4))
    
    return TensorMap(Tdag_data, space(T))
end