struct C4v <: AbstractSpaceGroup end

const C4v_σd_permutation = ((1, ), (3, 2, 5, 4))
const C4v_σv_permutation = ((1, ), (4, 3, 2, 5))
const C4v_R_permutation = ((1, ), (3, 4, 5, 2))

# in the order of σd, σv, R, see http://symmetry.jacobs-university.de/cgi-bin/group.cgi?group=404&option=4
const C4v_A1_reps = (1, 1, 1) 
const C4v_A2_reps = (-1, -1, 1) 
const C4v_B1_reps = (-1, 1, 1) 
const C4v_B2_reps = (1, -1, 1) 

function get_reps(::Type{C4v}, name::Symbol)
    (name == :A1) && return C4v_A1_reps
    (name == :A2) && return C4v_A2_reps
    (name == :B1) && return C4v_B1_reps
    (name == :B2) && return C4v_B2_reps
end
function get_perm(::Type{C4v}, name::Symbol)
    (name == :σd) && return C4v_σd_permutation
    (name == :σv) && return C4v_σv_permutation
    (name == :R) && return C4v_R_permutation
end