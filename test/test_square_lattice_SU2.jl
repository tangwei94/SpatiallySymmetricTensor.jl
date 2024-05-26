@testset "test T_1_3_A1" begin
    T1 = IPEPSC6v.T_1_3_A1()
    T2 = IPEPSC6v.T_1_3_A1_from_plain()
    @show λ = norm(T1) / norm(T2)

    @test norm(λ * T2 - T1) < 1e-12
end

@testset "test T_3_1_A1" begin
    T1 = IPEPSC6v.T_3_1_A1()
    T2 = IPEPSC6v.T_3_1_A1_from_plain()
    @show λ = norm(T1) / norm(T2)

    @test norm(λ * T2 - T1) < 1e-12
end