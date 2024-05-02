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

function mpo_hermicity(A)
    v_space = MPSKit._firstspace(A)
    function AA_transf(v)
        @tensor Tv[-1; -2] := A[-1 3; 4 1] * conj(A[-2 3; 4 2]) * v[1; 2]
        return Tv
    end
    function AĀ_transf(v)
        @tensor Tv[-1 -2] := A[-1 3; 4 1] * A[-2 4; 3 2] * v[1 2]
        return Tv
    end

    vaa = TensorMap(rand, ComplexF64, v_space, v_space)
    vaā = Tensor(rand, ComplexF64, v_space*v_space)

    aa, _ = eigsolve(AA_transf, vaa, 1, :LM) 
    aā, _ = eigsolve(AĀ_transf, vaā, 1, :LM) 

    aa = aa[1]
    aā = aā[1]

    return norm(aā / aa) 
end

function mpo_normality(A)
    v_space = MPSKit._firstspace(A)

    function ĀAĀA_transf(v)
        @tensor Tv[-1 -2; -3 -4] := conj(A[-3 7; 8 6]) * A[-1 7; 5 4] * conj(A[-4 3; 5 2]) * A[-2 3; 8 1] * v[1 4; 2 6]
        return Tv
    end
    function ĀAAĀ_transf(v)
        @tensor Tv[-1 -2; -3 -4] := conj(A[-3 7; 8 6]) * A[-1 7; 5 4] * A[-2 5; 3 2] * conj(A[-4 8; 3 1]) * v[4 2; 6 1]
        return Tv
    end
    
    v0 = TensorMap(rand, ComplexF64, v_space*v_space, v_space*v_space)

    āaāa, _ = eigsolve(ĀAĀA_transf, v0, 1, :LM) 
    āaaā, _ = eigsolve(ĀAAĀ_transf, v0, 1, :LM) 

    āaāa = āaāa[1]
    āaaā = āaaā[1]

    return norm(āaāa / āaaā) 
end