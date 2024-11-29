@testset "test T_1_3_A1" begin
    T1 = SpatiallySymmetricTensor.T_1_3_A1()
    T2 = SpatiallySymmetricTensor.T_1_3_A1_from_plain()
    _, ix = findmax(norm.(T1.data))
    λ = T1.data[ix] / T2.data[ix]

    @test norm(λ * T2 - T1) < 1e-12
end

@testset "test T_3_1_A1" begin
    T1 = SpatiallySymmetricTensor.T_3_1_A1()
    T2 = SpatiallySymmetricTensor.T_3_1_A1_from_plain()
    _, ix = findmax(norm.(T1.data))
    λ = T1.data[ix] / T2.data[ix]

    @test norm(λ * T2 - T1) < 1e-12
end

@testset "spin exchange" begin 
    Sleft, Sright = SpatiallySymmetricTensor.spin_exchange()
    @tensor SdotS[-1 -2; -3 -4] := Sleft[-1; -3 1] * Sright[1 -2; -4]
    SdotS_arr1 = reshape(convert(Array, SdotS), (4, 4))

    Sx = ComplexF64[0 0.5; 0.5 0]
    Sy = ComplexF64[0 -0.5im; 0.5im 0]
    Sz = ComplexF64[0.5 0; 0 -0.5]
    SdotS_arr2 = kron(Sx, Sx) + kron(Sy, Sy) + kron(Sz, Sz)

    @test norm(SdotS_arr1 - SdotS_arr2) < 1e-12
end

