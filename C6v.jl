# continue from C4v.jl

using Combinatorics

# fuse 6 spin (1/2)⊕(0) into (1/2), nocc = {3, 3} (long-range RVB, A_1^{2})
function _f(_is, _s1, _s2, _s3)
    _Szs = zeros(Union{Int, Rational{Int64}}, 6)
    _Szs[_is[1]] = _s1
    _Szs[_is[2]] = _s2
    _Szs[_is[3]] = _s3

    _Ss = zeros(Union{Int, Rational{Int64}}, 6)
    _Ss[_is] .= 1//2 
    return vstate(_Szs, _Ss) 
end
dat = vec(map(t -> _f(t...), Iterators.product(collect(combinations(1:6, 3)), [-1//2, 1//2], [-1//2, 1//2], [-1//2, 1//2])))
sort!(dat; by=Sz_order_score)

S2 = fS2(dat)
Sz = fSz(dat)

ΛS2, US2 = eigen(Hermitian(S2))
PS2 = US2[:, ΛS2 .≈ 0.75]
@show PS2' * S2 * PS2

ΛSz, USz = eigen(Hermitian(PS2' * Sz * PS2))
PSz = USz[:, ΛSz .≈ 1/2]

function σv_C6v(a::vstate)
    perm = [5, 4, 3, 2, 1, 6]
    return vstate(a.Szs[perm], a.Ss[perm])
end
function σd_C6v(a::vstate)
    perm = [6, 5, 4, 3, 2, 1]
    return vstate(a.Szs[perm], a.Ss[perm])
end

Rmat = f_space_op(dat, R)
ΛR, UR = eigen(PSz' * PS2' * Rmat * PS2 * PSz)
PR = UR[:, ΛR .≈ 1]
PRinv = inv(UR)[ΛR .≈ 1, :]

σv_mat = f_space_op(dat, σv_C6v)
Λσv, Uσv = eigen(PRinv * PSz' * PS2' * σv_mat * PS2 * PSz * PR)
Pσv = Uσv[:, Λσv .≈ 1]
Pσvinv = inv(Uσv)[Λσv .≈ 1, :]

σd_mat = f_space_op(dat, σd_C6v)
Pσvinv * PRinv * PSz' * PS2' * σd_mat * PS2 * PSz * PR * Pσv

sol = PS2 * PSz * PR * Pσv
sol1 = sol[:, 1] / norm(sol[:, 1])
sol2 = sol[:, 2] / norm(sol[:, 2])
sol3 = sol[:, 3] / norm(sol[:, 3])
soli = sol3;
for ix in eachindex(soli)
    if abs(soli[ix]) > 1e-12 && dat[ix].Ss[1:3] == [0//1, 0//1, 0//1] 
        printstyled("$(ix) $(sol[ix])\n"; color=:red, bold=true)
        printstyled("$(dat[ix])\n\n")
    end
end

σv_mat = f_space_op(dat, σv_C6v)
Λσv, Uσv = eigen(PRinv * PSz' * PS2' * σv_mat * PS2 * PSz * PR)
Pσv_A2 = Uσv[:, Λσv .≈ -1]
Pσvinv_A2 = inv(Uσv)[Λσv .≈ -1, :]

Pσvinv_A2 * PRinv * PSz' * PS2' * σd_mat * PS2 * PSz * PR * Pσv_A2

sol_A2 = PS2 * PSz * PR * Pσv_A2
sol1_A2 = sol_A2[:, 1] / norm(sol_A2[:, 1])
sol2_A2 = sol_A2[:, 2] / norm(sol_A2[:, 2])
sol3_A2 = sol_A2[:, 3] / norm(sol_A2[:, 3])
sol_i = sol1_A2;
for ix in eachindex(sol_i)
    if abs(sol_i[ix]) > 1e-12 && dat[ix].Ss[2:4] == [0//1, 0//1, 0//1] 
        printstyled("$(ix) $(sol_i[ix])\n"; color=:red, bold=true)
        printstyled("$(dat[ix])\n\n")
    end
end

# fuse 6 spin (1/2)⊕(0) into (1/2), nocc = {5, 1} (long-range RVB, A_1^{2})
function _f(_is, _s1, _s2, _s3, _s4, _s5)
    _Szs = zeros(Union{Int, Rational{Int64}}, 6)
    _Szs[_is[1]] = _s1
    _Szs[_is[2]] = _s2
    _Szs[_is[3]] = _s3
    _Szs[_is[4]] = _s4
    _Szs[_is[5]] = _s5

    _Ss = zeros(Union{Int, Rational{Int64}}, 6)
    _Ss[_is] .= 1//2 
    return vstate(_Szs, _Ss) 
end
dat = vec(map(t -> _f(t...), Iterators.product(collect(combinations(1:6, 5)), [-1//2, 1//2], [-1//2, 1//2], [-1//2, 1//2], [-1//2, 1//2], [-1//2, 1//2])))
sort!(dat; by=Sz_order_score)

S2 = fS2(dat)
Sz = fSz(dat)

ΛS2, US2 = eigen(Hermitian(S2))
PS2 = US2[:, ΛS2 .≈ 0.75]
@show PS2' * S2 * PS2

ΛSz, USz = eigen(Hermitian(PS2' * Sz * PS2))
PSz = USz[:, ΛSz .≈ 1/2]

Rmat = f_space_op(dat, R)
ΛR, UR = eigen(PSz' * PS2' * Rmat * PS2 * PSz)
PR = UR[:, ΛR .≈ 1]
PRinv = inv(UR)[ΛR .≈ 1, :]

σv_mat = f_space_op(dat, σv_C6v)
Λσv, Uσv = eigen(PRinv * PSz' * PS2' * σv_mat * PS2 * PSz * PR)
Pσv = Uσv[:, Λσv .≈ 1]
Pσvinv = inv(Uσv)[Λσv .≈ 1, :]

σd_mat = f_space_op(dat, σd_C6v)
Pσvinv * PRinv * PSz' * PS2' * σd_mat * PS2 * PSz * PR * Pσv

sol = PS2 * PSz * PR * Pσv
sol1 = sol[:, 1] / norm(sol[:, 1])
sol2 = sol[:, 2] / norm(sol[:, 2])
sol3 = sol[:, 3] / norm(sol[:, 3])
soli = sol3;
for ix in eachindex(soli)
    if abs(soli[ix]) > 1e-12 
        printstyled("$(ix) $(sol[ix])\n"; color=:red, bold=true)
        printstyled("$(dat[ix])\n\n")
    end
end

σv_mat = f_space_op(dat, σv_C6v)
Λσv, Uσv = eigen(PRinv * PSz' * PS2' * σv_mat * PS2 * PSz * PR)
Pσv_A2 = Uσv[:, Λσv .≈ -1]
Pσvinv_A2 = inv(Uσv)[Λσv .≈ -1, :]

Pσvinv_A2 * PRinv * PSz' * PS2' * σd_mat * PS2 * PSz * PR * Pσv_A2

sol_A2 = PS2 * PSz * PR * Pσv_A2
sol1_A2 = sol_A2[:, 1] / norm(sol_A2[:, 1])
sol2_A2 = sol_A2[:, 2] / norm(sol_A2[:, 2])
sol_i = sol1_A2;
for ix in eachindex(sol_i)
    if abs(sol_i[ix]) > 1e-12 
        printstyled("$(ix) $(sol_i[ix])\n"; color=:red, bold=true)
        printstyled("$(dat[ix])\n\n")
    end
end

