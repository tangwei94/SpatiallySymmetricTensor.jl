struct C4 <: AbstractPointGroup end

# in the order of σd, σv, R; see http://symmetry.jacobs-university.de/cgi-bin/group.cgi?group=204&option=4
# ignore representationq of σd and σv for C4
const C4_A_reps = (0, 0, 1) 
const C4_B_reps = (0, 0, -1) 

const C4_R1 = ((1, ), (3, 4, 5, 2))
const C4_R2 = ((1, ), (5, 2, 3, 4))

function get_reps(::C4, name::Symbol)
    (name == :A) && return C4_A_reps
    (name == :B) && return C4_B_reps
end
function get_perm(::C4, name::Symbol)
    (name == :σd) && return [] # empty diagonal reflections
    (name == :σv) && return [] # empty horizontal and vertical reflections
    (name == :R) && return [C4_R1, C4_R2]
end
