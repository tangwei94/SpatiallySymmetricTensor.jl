using LinearAlgebra 

struct vstate
    Szs::Vector{Union{Int64, Rational{Int64}}}
    Ss::Vector{Union{Int64, Rational{Int64}}}
end

Sz_order_score(a::vstate) = sum(10 .^ (1:length(a.Szs)) .* a.Szs)

# only consider the case of V = 1//2 ⊕ 0
function fS2(_dat::Vector{vstate})
    N = length(_dat)
    S2 = zeros(N, N)
    for (_i, _si) in enumerate(_dat), (_j, _sj) in enumerate(_dat)
        _si.Ss == _sj.Ss &&  
        sum(_si.Szs) == sum(_sj.Szs) &&  
        sum(abs.(_si.Szs .- _sj.Szs)) == 2 && 
        (S2[_i, _j] = 1)
    end
    for (_i, _si) in enumerate(_dat)
        S2[_i, _i] = 2*sum(_si.Ss)/2 + sum(_si.Szs)^2 
    end
    return S2
end

function fSz(_dat::Vector{vstate})
    N = length(_dat)
    Sz = zeros(N, N)
    for (_i, _si) in enumerate(_dat)
        Sz[_i, _i] = sum(_si.Szs) 
    end
    return Sz
end

## example 1: fuse 3 spin 1/2
function _f(_s1, _s2, _s3)  
    return vstate([_s1, _s2, _s3], [1//2, 1//2, 1//2])
end
dat = vec(map(t -> _f(t...), Iterators.product([-1//2, 1//2], [-1//2, 1//2], [-1//2, 1//2])))
sort!(dat; by=Sz_order_score)

S2 = fS2(dat)
Sz = fSz(dat)

ΛS2, US2 = eigen(Hermitian(S2))
P1 = US2[:, 1:4]
@show P1' * S2 * P1

ΛSz, USz = eigen(Hermitian(P1' * Sz * P1))
P1 * USz[:, 1]
dat[2]
dat[3]

P1 * USz[:, 2]
dat[2]
dat[3]
dat[5]

## example 2: fuse 4 spin (1/2)⊕(0) into 1/2, nocc = {3, 1} (long-range RVB, A_1^{2})
function _f(_i, _s1, _s2, _s3) # _i: location of spin 0, 
    _Szs = zeros(Union{Int, Rational{Int64}}, 4)
    _Szs[_i] = 0
    _Szs[(_i + 0) % 4 + 1] = _s1
    _Szs[(_i + 1) % 4 + 1] = _s2
    _Szs[(_i + 2) % 4 + 1] = _s3

    _Ss = fill(1//2, 4)
    _Ss[_i] = 0
    return vstate(_Szs, _Ss)
end
dat = vec(map(t -> _f(t...), Iterators.product(1:4, [-1//2, 1//2], [-1//2, 1//2], [-1//2, 1//2])))
sort!(dat; by=Sz_order_score)

S2 = fS2(dat)
Sz = fSz(dat)

ΛS2, US2 = eigen(Hermitian(S2))
PS2 = US2[:, ΛS2 .≈ 0.75]
@show PS2' * S2 * PS2

ΛSz, USz = eigen(Hermitian(PS2' * Sz * PS2))
PSz = USz[:, ΛSz .≈ 1/2]

function R(a::vstate)
    return vstate([a.Szs[end]; a.Szs[1:end-1]], [a.Ss[end]; a.Ss[1:end-1]])
end
function σv_C4v(a::vstate)
    perm = [3, 2, 1, 4]
    return vstate(a.Szs[perm], a.Ss[perm])
end
function σd_C4v(a::vstate)
    perm = [4, 3, 2, 1]
    return vstate(a.Szs[perm], a.Ss[perm])
end

function f_space_op(_dat::Vector{vstate}, _op::Function)
    @assert issorted(_dat; by=Sz_order_score)
    N = length(_dat)
    mat = zeros(N, N)
    for (_i, _si) in enumerate(_dat)
        _j = searchsortedfirst(_dat, _op(_si); by=Sz_order_score)
        mat[_j, _i] = 1
    end
    return mat
end

Rmat = f_space_op(dat, R)
ΛR, UR = eigen(PSz' * PS2' * Rmat * PS2 * PSz)
PR = UR[:, ΛR .≈ 1]
PRinv = inv(UR)[ΛR .≈ 1, :]

σv_mat = f_space_op(dat, σv_C4v)
Λσv, Uσv = eigen(PRinv * PSz' * PS2' * σv_mat * PS2 * PSz * PR)
Pσv = Uσv[:, Λσv .≈ 1]
Pσvinv = inv(Uσv)[Λσv .≈ 1, :]

σd_mat = f_space_op(dat, σd_C4v)
Pσvinv * PRinv * PSz' * PS2' * σd_mat * PS2 * PSz * PR * Pσv

sol = PS2 * PSz * PR * Pσv
sol /= norm(sol)
for ix in eachindex(sol)
    if abs(sol[ix]) > 1e-12 
        printstyled("$(ix) $(sol[ix])\n"; color=:red, bold=true)
        printstyled("$(dat[ix])\n\n")
    end
end

σv_mat = f_space_op(dat, σv_C4v)
Λσv, Uσv = eigen(PRinv * PSz' * PS2' * σv_mat * PS2 * PSz * PR)
Pσv_A2 = Uσv[:, Λσv .≈ -1]
Pσvinv_A2 = inv(Uσv)[Λσv .≈ -1, :]

Pσvinv_A2 * PRinv * PSz' * PS2' * σd_mat * PS2 * PSz * PR * Pσv_A2

sol_A2 = PS2 * PSz * PR * Pσv_A2
sol_A2 /= norm(sol_A2)
for ix in eachindex(sol_A2)
    if abs(sol[ix]) > 1e-12 
        printstyled("$(ix) $(sol[ix])\n"; color=:red, bold=true)
        printstyled("$(dat[ix])\n\n")
    end
end
