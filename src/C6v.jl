struct C6v <: AbstractPointGroup end

const C6v_R1 = ((1, ), (3, 4, 5, 6, 7, 2))
const C6v_R2 = ((1, ), (7, 2, 3, 4, 5, 6))
const C6v_σd1 = ((1, ), (7, 6, 5, 4, 3, 2))
const C6v_σd2 = ((1, ), (3, 2, 7, 6, 5, 4))
const C6v_σd3 = ((1, ), (5, 4, 3, 2, 7, 6))
const C6v_σv1 = ((1, ), (6, 5, 4, 3, 2, 7))
const C6v_σv2 = ((1, ), (2, 7, 6, 5, 4, 3))
const C6v_σv3 = ((1, ), (4, 3, 2, 7, 6, 5))

# in the order of σd, σv, R, see http://symmetry.jacobs-university.de/cgi-bin/group.cgi?group=406&option=4
const C6v_A1_reps = (1, 1, 1) 
const C6v_A2_reps = (-1, -1, 1) 
const C6v_B1_reps = (-1, 1, -1) 
const C6v_B2_reps = (1, -1, -1) 

function get_reps(::C6v, name::Symbol)
    (name == :A1) && return C6v_A1_reps
    (name == :A2) && return C6v_A2_reps
    (name == :B1) && return C6v_B1_reps
    (name == :B2) && return C6v_B2_reps
end
function get_perm(::C6v, name::Symbol)
    (name == :σd) && return [C6v_σd1, C6v_σd2, C6v_σd3] 
    (name == :σv) && return [C6v_σv1, C6v_σv2, C6v_σv3]
    (name == :R) && return [C6v_R1, C6v_R2]
end
