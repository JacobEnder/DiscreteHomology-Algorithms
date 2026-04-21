#===================================================================#
# Embeddings.jl
# Injective graph homomorphisms (subgraph monomorphisms) from a small
# source graph into a target graph.
#
# The source graphs here are tiny (<= 8 vertices), so a backtracking
# search with adjacency pruning is adequate.
#===================================================================#

"""Adjacency sets (no loops) for fast intersection."""
function adjacency_sets_no_loops(G::SimpleGraph)
    adj = adjacency_no_loops(G)
    return Dict(v => Set(nbs) for (v, nbs) in adj)
end

"""Degree (excluding loops) of each vertex."""
function degrees_no_loops(G::SimpleGraph)
    adj = adjacency_no_loops(G)
    return [length(adj[v]) for v in G.vertices]
end

"""
    find_injective_homs(src, tgt)

Return a vector of injective graph homomorphisms src -> tgt.

Each map is returned as a Vector{Int} of length |V(src)|, where
map[v] is the image of src vertex v in tgt.
"""
function find_injective_homs(src::SimpleGraph, tgt::SimpleGraph)
    n = length(src.vertices)
    m = length(tgt.vertices)
    n == 0 && return Vector{Vector{Int}}([Int[]])
    n > m && return Vector{Vector{Int}}()  # cannot inject

    src_adj = adjacency_sets_no_loops(src)
    tgt_adj = adjacency_sets_no_loops(tgt)
    src_deg = degrees_no_loops(src)
    tgt_deg = degrees_no_loops(tgt)

    # Order src vertices by descending degree for pruning
    order = sortperm(src_deg; rev=true)

    pos_of = zeros(Int, n)
    for (p, v) in enumerate(order)
        pos_of[v] = p
    end

    # For each position in 'order', record which earlier vertices are neighbors
    earlier_neighbors = Vector{Vector{Int}}(undef, n)
    for pos in 1:n
        v = order[pos]
        nbs = Int[]
        for prevpos in 1:(pos-1)
            u = order[prevpos]
            (u in src_adj[v]) && push!(nbs, u)
        end
        earlier_neighbors[pos] = nbs
    end

    # Candidate targets by degree: image vertex must have degree >= preimage degree
    cand_by_deg = [Int[] for _ in 1:n]
    for pos in 1:n
        v = order[pos]
        dv = src_deg[v]
        cand_by_deg[pos] = [w for w in tgt.vertices if tgt_deg[w] >= dv]
    end

    current = fill(0, n)   # mapping in ordered coordinates: current[pos] = tgt vertex
    used = falses(m)
    maps = Vector{Vector{Int}}()

    function backtrack(pos::Int)
        if pos > n
            # unpermute into original vertex order (1..n)
            f = zeros(Int, n)
            for p in 1:n
                f[order[p]] = current[p]
            end
            push!(maps, f)
            return
        end

        # adjacency constraint candidates
        candidates = cand_by_deg[pos]
        if !isempty(earlier_neighbors[pos])
            # intersect neighbors of already assigned images
            # start with all candidates; refine by each constraint
            filtered = Set(candidates)
            for u in earlier_neighbors[pos]
                img_u = current[pos_of[u]]
                filtered = intersect(filtered, tgt_adj[img_u])
                isempty(filtered) && return
            end
            candidates = sort!(collect(filtered))
        end

        for w in candidates
            used[w] && continue
            # Check adjacency to earlier neighbors explicitly (since candidates were intersected,
            # but we also need to ensure non-neighbors do not force edges; for homomorphisms we
            # only need to preserve edges, not non-edges.)
            ok = true
            for u in earlier_neighbors[pos]
                img_u = current[pos_of[u]]
                (w in tgt_adj[img_u]) || (w == img_u) || (ok=false)
                ok || break
            end
            ok || continue

            current[pos] = w
            used[w] = true
            backtrack(pos + 1)
            used[w] = false
            current[pos] = 0
        end
    end

    backtrack(1)
    return maps
end
