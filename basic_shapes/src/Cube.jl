#===================================================================#
# Cube.jl
# Discrete 3-cube Q^3 utilities.
#
# Vertex labeling (1..8) follows binary triples:
#   1=(0,0,0) 2=(0,0,1) 3=(0,1,0) 4=(0,1,1)
#   5=(1,0,0) 6=(1,0,1) 7=(1,1,0) 8=(1,1,1)
#
# Faces are ordered 4-cycles (a,b,c,d) tracing the square boundary.
#===================================================================#

"""Binary coordinates of Q^3 vertices (1-indexed)."""
const Q3_COORDS = NTuple{3,Int}[
    (0, 0, 0),
    (0, 0, 1),
    (0, 1, 0),
    (0, 1, 1),
    (1, 0, 0),
    (1, 0, 1),
    (1, 1, 0),
    (1, 1, 1),
]

"""Map a binary triple to vertex index 1..8."""
const Q3_INDEX = Dict{NTuple{3,Int},Int}(Q3_COORDS[i] => i for i in 1:8)

"""Flip the d-th coordinate (d in 1..3) of vertex i."""
function q3_flip(i::Int, d::Int)
    c = Q3_COORDS[i]
    newc = d == 1 ? (1 - c[1], c[2], c[3]) :
           d == 2 ? (c[1], 1 - c[2], c[3]) :
           d == 3 ? (c[1], c[2], 1 - c[3]) :
           error("dimension d must be 1..3")
    return Q3_INDEX[newc]
end

"""Precomputed flip table FLIP[d][i] = i with d-th bit flipped."""
const Q3_FLIP = [ [q3_flip(i, d) for i in 1:8] for d in 1:3 ]

"""Build the reflexive graph Q^3 as a SimpleGraph."""
function build_q3()
    vertices = collect(1:8)
    raw_edges = Tuple{Int,Int}[]
    for i in 1:8
        for j in (i+1):8
            ci, cj = Q3_COORDS[i], Q3_COORDS[j]
            diff = (ci[1] != cj[1]) + (ci[2] != cj[2]) + (ci[3] != cj[3])
            diff == 1 && push!(raw_edges, (i, j))
        end
    end
    edges = make_symmetric_reflexive(raw_edges, vertices)
    adjacency = build_adjacency(vertices, edges)
    return SimpleGraph("Q3", vertices, edges, adjacency)
end

"""The 6 faces of Q^3 as ordered square boundaries."""
const Q3_FACES = NTuple{4,Int}[
    (1, 2, 4, 3),   # x1=0
    (5, 6, 8, 7),   # x1=1
    (1, 2, 6, 5),   # x2=0
    (3, 4, 8, 7),   # x2=1
    (1, 3, 7, 5),   # x3=0
    (2, 4, 8, 6),   # x3=1
]

"""Normalize a 3-cycle to an ordered tuple (i<j<k)."""
normalize_3_cycle(a::Int, b::Int, c::Int) = Tuple(sort((a, b, c)))

"""Normalize a 4-cycle to canonical form (min over rotations/reflections)."""
function normalize_4_cycle(a::Int, b::Int, c::Int, d::Int)
    reps = NTuple{4,Int}[
        (a, b, c, d), (b, c, d, a), (c, d, a, b), (d, a, b, c),
        (a, d, c, b), (d, c, b, a), (c, b, a, d), (b, a, d, c),
    ]
    return minimum(reps)
end

"""Normalize either a 3-tuple or 4-tuple cycle."""
normalize_cycle(c::NTuple{3,Int}) = normalize_3_cycle(c...)
normalize_cycle(c::NTuple{4,Int}) = normalize_4_cycle(c...)

"""
    face_cycle_under_map(face, f)

Push a Q^3 face (ordered 4-cycle) through a vertex map f.

`f` is a vector indexed by Q^3 vertices; `f[v]` is the image of vertex v.

We collapse consecutive duplicates cyclically; if the resulting cycle has
3 or 4 distinct vertices we return a tuple. Otherwise return `nothing`.
"""
function face_cycle_under_map(face::NTuple{4,Int}, f::AbstractVector{Int})
    mapped = [f[v] for v in face]

    changed = true
    while changed
        changed = false
        n = length(mapped)
        n == 0 && return nothing
        new_seq = Int[]
        for i in 1:n
            j = (i % n) + 1
            if mapped[i] != mapped[j]
                push!(new_seq, mapped[i])
            else
                changed = true
            end
        end
        mapped = new_seq
    end

    length(mapped) < 3 && return nothing
    length(unique(mapped)) != length(mapped) && return nothing

    if length(mapped) == 3
        return (mapped[1], mapped[2], mapped[3])
    elseif length(mapped) == 4
        return (mapped[1], mapped[2], mapped[3], mapped[4])
    end
    return nothing
end

"""
    is_degenerate_3cube_map(f)

BCW degeneracy for a 3-cube map: the map is constant in some coordinate.
Equivalently, for some dimension d, f(i) == f(flip_d(i)) for all i.
"""
function is_degenerate_3cube_map(f::AbstractVector{Int})
    @assert length(f) == 8
    for d in 1:3
        ok = true
        for i in 1:8
            if f[i] != f[Q3_FLIP[d][i]]
                ok = false
                break
            end
        end
        ok && return true
    end
    return false
end
