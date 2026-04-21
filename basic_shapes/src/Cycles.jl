#===================================================================#
# Cycles.jl
# Enumerate 3-cycles and 4-cycles in a simple undirected graph
# (loops are ignored).
#
# 3-cycles are returned as sorted tuples (i<j<k).
# 4-cycles are returned in the canonical form used by normalize_4_cycle.
#===================================================================#

using SparseArrays

"""Find all triangles in the graph given a no-loop adjacency dict."""
function detect_3_cycles(adj::Dict{Int, Vector{Int}}, n_vertices::Int)
    I_arr, J_arr = Int[], Int[]
    for (v, nbs) in adj
        for u in nbs
            push!(I_arr, v)
            push!(J_arr, u)
        end
    end
    isempty(I_arr) && return NTuple{3,Int}[]
    A = sparse(I_arr, J_arr, ones(Int, length(I_arr)), n_vertices, n_vertices)

    tris = Set{NTuple{3,Int}}()
    for i in 1:n_vertices
        haskey(adj, i) || continue
        nbs = adj[i]
        for (idx, j) in enumerate(nbs)
            j <= i && continue
            for k in nbs[(idx+1):end]
                (k > j && A[j, k] > 0) || continue
                push!(tris, (i, j, k))
            end
        end
    end
    return sort!(collect(tris))
end

"""Find all 4-cycles using the common-neighbors method."""
function detect_4_cycles(adj::Dict{Int, Vector{Int}}, n_vertices::Int)
    adj_sets = Dict(v => Set(nbs) for (v, nbs) in adj)
    cycles = Set{NTuple{4,Int}}()
    for u in 1:n_vertices
        haskey(adj_sets, u) || continue
        for w in (u+1):n_vertices
            haskey(adj_sets, w) || continue
            common = collect(intersect(adj_sets[u], adj_sets[w]))
            for i in 1:length(common)
                for j in (i+1):length(common)
                    c1, c2 = common[i], common[j]
                    push!(cycles, normalize_4_cycle(u, c1, w, c2))
                end
            end
        end
    end
    return sort!(collect(cycles))
end

"""Build lookup dict cycle -> C_2 index (triangles first, then squares)."""
function build_cycle_to_idx(triangles::Vector{NTuple{3,Int}}, four_cycles::Vector{NTuple{4,Int}})
    dict = Dict{Any, Int}()
    for (i, tri) in enumerate(triangles)
        dict[tri] = i
    end
    offset = length(triangles)
    for (i, sq) in enumerate(four_cycles)
        dict[sq] = offset + i
    end
    return dict
end
