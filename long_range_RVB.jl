# https://doi.org/10.1103/PhysRevLett.111.037202
# https://doi.org/10.1103/PhysRevB.94.205124

using TensorKit, LinearAlgebra, MPSKit, KrylovKit
using JLD2
using Revise
using IPEPSC6v

V = SU2Space(1//2=>1, 0=>1)
P = SU2Space(1//2=>1)
T = TensorMap(zeros, ComplexF64, P, V^4)

# a projector to the subspace with (1//2 ⊕ 0 ⊕ 0 ⊕ 0 -> 1//2) (short-range RVB)
P_nocc_1_3 = begin
    _condition(f1, f2) = length(findall(rep-> rep == SU2Irrep(1//2), f2.uncoupled)) == 1
    selector(T, _condition)
end

# a projector to the subspace with (1//2 ⊕ 1//2 ⊕ 1//2 ⊕ 0 -> 1//2) (long-range RVB)
P_nocc_3_1 = begin
    _condition(f1, f2) = length(findall(rep-> rep == SU2Irrep(1//2), f2.uncoupled)) == 3
    selector(T, _condition)
end

# matrices for spatial operations
R_mat = spatial_operation(T, ((1, ), (3, 4, 5, 2)))
σd_mat = spatial_operation(T, ((1, ), (3, 2, 5, 4)))
σv_mat = spatial_operation(T, ((1, ), (4, 3, 2, 5)))

##### nocc= {1, 3}, TABLE VII in PRB 94, 205124 (2016)
T_1_3_A1 = begin
    Λσv_1_3, Uσv_1_3 = eigen(Hermitian(P_nocc_1_3' * σv_mat * P_nocc_1_3))
    Pσv_1_3 = Uσv_1_3[:, Λσv_1_3 .≈ 1]

    Λσd_1_3, Uσd_1_3 = eigen(Hermitian(Pσv_1_3' * P_nocc_1_3' * σd_mat * P_nocc_1_3 * Pσv_1_3))
    Pσd_1_3 = Uσd_1_3[:, Λσd_1_3 .≈ 1]

    ΛR_1_3, UR_1_3 = eigen(Pσd_1_3' * Pσv_1_3' * P_nocc_1_3' * R_mat * P_nocc_1_3 * Pσv_1_3 * Pσd_1_3)
    PR_1_3 = UR_1_3[:, ΛR_1_3 .≈ 1]

    sol_1_3 = vec(P_nocc_1_3 * Pσv_1_3 * Pσd_1_3 * PR_1_3)
    T_1_3_A1 = set_data_by_vector(T, sol_1_3)
end
##### nocc= {3, 1}, TABLE IX.  in PRB 94, 205124 (2016)
T_3_1_A1 = begin
    Λσv_3_1, Uσv_3_1 = eigen(Hermitian(P_nocc_3_1' * σv_mat * P_nocc_3_1))
    Pσv_3_1 = Uσv_3_1[:, Λσv_3_1 .≈ 1]

    Λσd_3_1, Uσd_3_1 = eigen(Hermitian(Pσv_3_1' * P_nocc_3_1' * σd_mat * P_nocc_3_1 * Pσv_3_1))
    Pσd_3_1 = Uσd_3_1[:, Λσd_3_1 .≈ 1]

    ΛR_3_1, UR_3_1 = eigen(Pσd_3_1' * Pσv_3_1' * P_nocc_3_1' * R_mat * P_nocc_3_1 * Pσv_3_1 * Pσd_3_1)
    PR_3_1 = UR_3_1[:, ΛR_3_1 .≈ 1]

    sol_3_1 = vec(P_nocc_3_1 * Pσv_3_1 * Pσd_3_1 * PR_3_1)
    T_3_1_A1 = set_data_by_vector(T, sol_3_1)
end

# get the two symmetric tensors
T_1_3_A1 = T_1_3_A1 / norm(T_1_3_A1)
T_3_1_A1 = T_3_1_A1 / norm(T_3_1_A1)

λ = parse(Float64, ARGS[1])
A = T_1_3_A1 + λ * T_3_1_A1;
B = Tensor(zeros, ComplexF64, V*V);
B.data.values[1] .= [1.0, sqrt(2)] ;
Bdata = convert(Array, B);

δ = isomorphism(fuse(V'*V), V'*V);
δ = permute(δ, (1, 2), (3, ));
@tensor TA[-1 -2; -3 -4] := δ[-1 2; 1] * δ[-2 4; 3] * conj(δ[-3 5; 6]) * conj(δ[-4 7; 8]) * A[9; 4 2 6 8] * conj(A[9; 3 1 5 7]);
@tensor TB[-1; -2] := δ[-1 1; 2] * conj(δ[-2 4; 3]) * B[2 4] * conj(B[1 3]);

@tensor Tfull[-1 -2; -3 -4] := TA[-1 -2; 1 2] * TB[1; -3] * TB[2; -4];
IPEPSC6v.mpo_hermicity(Tfull)
IPEPSC6v.mpo_normality(Tfull)

#U1, _, _ = tsvd(Tfull, (1, ), (2, 3, 4); trunc=truncerr(1e-12)); 
#U2, _, _ = tsvd(Tfull, (2, ), (1, 3, 4); trunc=truncerr(1e-12));
#@tensor Tfull1[-1 -2; -3 -4] := Tfull[1 2; 3 4] * U1[4; -4] * U1'[-1; 1] * U2[3; -3] * U2'[-2; 2]; 
#
#IPEPSC6v.mpo_hermicity(Tfull1)
#IPEPSC6v.mpo_normality(Tfull1)

ψi = InfiniteMPS([fuse(V'*V)], [fuse(V'*V)])

𝕋full = DenseMPO([Tfull]) 
let ψ1 = ψi, ψ2 = ψi
    ψ1 = ψi
    for ix in 1:20 
        ψ1 = changebonds(𝕋full * ψ1, SvdCut(truncdim(100))) 
        @show ix, domain(ψ1.CR[1]) 
    end 
    ψ2, _, _ = leading_boundary(ψ1, 𝕋full, VUMPS(tol_galerkin=1e-12, maxiter=1000)); 
    @save "data/long_range_RVB_lambda$(λ)_chi100.jld2" ψ1 ψ2
end

@load "data/long_range_RVB_lambda$(λ)_chi100.jld2" ψ1 ψ2

# === compute ground state energy

# spin exchange
VA = SU2Space(1 => 1)
Sleft = TensorMap(ones, ComplexF64, P, P*VA) * sqrt(3/4)
Sright = -TensorMap(ones, ComplexF64, VA*P, P) * sqrt(3/4)

@tensor TSl[-1 -2; -3 -4 -5] := δ[-1 2; 1] * δ[-2 4; 3] * conj(δ[-3 5; 6]) * conj(δ[-4 7; 8]) * A[10; 4 2 6 8] * conj(A[9; 3 1 5 7]) * Sleft[9; 10 -5];
@tensor TSr[-1 -2 -3; -4 -5] := δ[-2 2; 1] * δ[-3 4; 3] * conj(δ[-4 5; 6]) * conj(δ[-5 7; 8]) * A[10; 4 2 6 8] * conj(A[9; 3 1 5 7]) * Sright[-1 9; 10];
@tensor TSleft[-1 -2; -3 -4 -5] := TSl[-1 -2; 1 2 -5] * TB[1; -3] * TB[2; -4];
@tensor TSright[-1 -2; -3 -4 -5] := TSr[-1 -2 -3; 1 2] * TB[1; -4] * TB[2; -5];

# transfer matrix
ψA = ψ2.AL[1]
function transfer_R(v)
    @tensor Tv[-1 -2 -3; -4] := ψA[-3 2; 1] * Tfull[-2 4; 2 3] * Tfull[-1 7; 4 5] * conj(ψA[-4 7; 6]) * v[5 3 1; 6]
    return Tv
end
function transfer_L(v)
    @tensor Tv[-1; -2 -3 -4] := ψA[1 2; -4] * Tfull[3 4; 2 -3] * Tfull[5 7; 4 -2] * conj(ψA[6 7; -1]) * v[6; 5 3 1]
    return Tv
end

Vψ = MPSKit._firstspace(ψA)
VT = MPSKit._firstspace(Tfull)

vr0 = TensorMap(rand, ComplexF64, VT*VT*Vψ, Vψ)
vl0 = TensorMap(rand, ComplexF64, Vψ, VT*VT*Vψ)

λr, Er = eigsolve(transfer_R, vr0, 1, :LM);
λr = λr[1];
Er = Er[1];

λl, El = eigsolve(transfer_L, vl0, 1, :LM);
λl = λl[1];
El = El[1];

Elr_ovlp = tr(El * Er)

@tensor E1 = El[2; 9 5 1] * ψA[1 6; 3] * TSleft[5 10; 6 7 12] * TSright[12 9 8; 10 11] * conj(ψA[2 8; 4]) * Er[11 7 3; 4] ;
E1 = E1 / λl / Elr_ovlp

@tensor E2 = El[2; 6 4 1] * ψA[1 3; 8] * TSleft[4 5; 3 11 13] * Tfull[6 7; 5 14] * conj(ψA[2 7; 9]) * ψA[8 10; 16] * TSright[13 11 12; 10 18] * Tfull[14 15; 12 19] * conj(ψA[9 15; 17]) * Er[19 18 16; 17];
E2 = E2 / λl^2 / Elr_ovlp

@tensor E3 = El[2; 6 4 1] * ψA[1 3; 8] * TSleft[4 5; 3 11 13] * Tfull[6 7; 5 14] * conj(ψA[2 7; 9]) * ψA[8 10; 16] * Tfull[11 12; 10 18] * TSright[13 14 15; 12 19] * conj(ψA[9 15; 17]) * Er[19 18 16; 17];
E3 = E3 / λl^2 / Elr_ovlp

@tensor E4 = El[2; 6 4 1] * ψA[1 3; 8] * TSleft[6 7; 5 14 13] * Tfull[4 5; 3 11] * conj(ψA[2 7; 9]) * ψA[8 10; 16] * Tfull[14 15; 12 19] * TSright[13 11 12; 10 18] * conj(ψA[9 15; 17]) * Er[19 18 16; 17];
E4 = E4 / λl^2 / Elr_ovlp
 
J2 = 0.5
@show E1, E2, E3, E4
E1, E2, E3, E4 = real(E1), real(E2), real(E3), real(E4)
Etot = E1 + E2 + J2*(E3 + E4)

io = open("tmpdata.txt", "a");
write(io, "$(λ) $(E1) $(E2) $(E3) $(E4) $(Etot)\n")
close(io)