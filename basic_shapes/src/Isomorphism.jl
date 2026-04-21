#===================================================================#
# Isomorphism.jl
# Brute-force graph isomorphism for very small graphs (<=8 vertices).
#
# We use degree-sequence pruning + backtracking on permutations.
# This is not intended as a general-purpose GI solver.
#===================================================================#

"""Boolean adjacency matrix without self-loops."""
function adjacency_matrix_bool(G::SimpleGraph)
    n = length(G.vertices)
    A = falses(n, n)
    for (u, v) in G.edges
        u == v && continue
        A[u, v] = true
    end
    return A
end

"""Sorted degree sequence excluding loops."""
function degree_sequence(G::SimpleGraph)
    degs = [count(u -> u != v, G.adjacency[v]) for v in G.vertices]
    sort!(degs)
    return degs
end

"""
    are_isomorphic(G1, G2)

Return (true, perm) if G1 and G2 are isomorphic, where perm maps
vertex i of G1 to vertex perm[i] of G2.
"""
function are_isomorphic(G1::SimpleGraph, G2::SimpleGraph)
    n1, n2 = length(G1.vertices), length(G2.vertices)
    n1 == n2 || return (false, Int[])
    n = n1

    degree_sequence(G1) == degree_sequence(G2) || return (false, Int[])

    A1 = adjacency_matrix_bool(G1)
    A2 = adjacency_matrix_bool(G2)

    deg1 = [count(u -> u != v, G1.adjacency[v]) for v in G1.vertices]
    deg2 = [count(u -> u != v, G2.adjacency[v]) for v in G2.vertices]

    g1_order = sortperm(deg1)  # assign low-degree first

    deg_to_g2 = Dict{Int, Vector{Int}}()
    for v in 1:n
        d = deg2[v]
        push!(get!(deg_to_g2, d, Int[]), v)
    end

    perm = zeros(Int, n)
    used = falses(n)

    function backtrack(depth::Int)
        depth > n && return true
        v = g1_order[depth]
        d = deg1[v]
        candidates = get(deg_to_g2, d, Int[])
        isempty(candidates) && return false
        for w in candidates
            used[w] && continue

            # adjacency consistency with already assigned vertices
            ok = true
            for prev in 1:(depth-1)
                u = g1_order[prev]
                pu = perm[u]
                if A1[u, v] != A2[pu, w]
                    ok = false
                    break
                end
            end
            ok || continue

            perm[v] = w
            used[w] = true
            backtrack(depth + 1) && return true
            perm[v] = 0
            used[w] = false
        end
        return false
    end

    found = backtrack(1)
    return (found, found ? copy(perm) : Int[])
end

"""Apply a vertex permutation perm (v -> perm[v]) to a 3- or 4-cycle."""
function permute_cycle(cycle::NTuple{3,Int}, perm::Vector{Int})
    return normalize_cycle((perm[cycle[1]], perm[cycle[2]], perm[cycle[3]]))
end
function permute_cycle(cycle::NTuple{4,Int}, perm::Vector{Int})
    return normalize_cycle((perm[cycle[1]], perm[cycle[2]], perm[cycle[3]], perm[cycle[4]]))
end
