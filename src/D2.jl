struct D2 <: AbstractPointGroup end

# in the order of σd, σv, R; see http://symmetry.jacobs-university.de/cgi-bin/group.cgi?group=302&option=4
# ignore representation of σd for D2
const D2_A_reps = (0, 1, 1) 
const D2_B1_reps = (0, -1, 1)

const D2_σv1 = ((1, ), (4, 3, 2, 5))
const D2_σv2 = ((1, ), (2, 5, 4, 3))

function get_reps(::D2, name::Symbol)
    (name == :A) && return D2_A_reps
    (name == :B1) && return D2_B1_reps
end
function get_perm(::D2, name::Symbol)
    (name == :σd) && return [] # empty diagonal reflections
    (name == :σv) && return [D2_σv1, D2_σv2]
    (name == :R) && return [] # empty rotations
end
