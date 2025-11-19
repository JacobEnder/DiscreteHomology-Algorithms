# Functions to detect and report all simple 3-cycles and 4-cycles in an undirected graph. 
# We use the natural ordering of vertices (assuming they are comparable)
# to ensure each triangle is reported exactly once.

module CycleDetection

using DataStructures
using SparseArrays
using Base.Threads

export detect_3_cycles, detect_4_cycles

function detect_3_cycles(adj)

    vertices = sort(collect(keys(adj)))
    n = length(vertices)
    vertex_to_idx = Dict(v => i for (i, v) in enumerate(vertices))
    
    # Build sparse adjacency matrix
    I = Int[]
    J = Int[]
    
    for v in vertices
        i = vertex_to_idx[v]
        for u in adj[v]
            if haskey(vertex_to_idx, u)
                j = vertex_to_idx[u]
                push!(I, i)
                push!(J, j)
            end
        end
    end
    
    A = sparse(I, J, ones(Int, length(I)), n, n)

    # Pre-compute neighbor lists for all vertices 
    neighbors_list = Vector{Vector{Int}}(undef, n)
    for i in 1:n
        neighbors_list[i] = findall(x -> x > 0, A[i, :])
    end

    # Thread-safe triangle collection
    triangles = Vector{Vector{Tuple{Any,Any,Any}}}(undef, nthreads())
    for i in 1:nthreads()
        triangles[i] = []
    end

    @threads for i in 1:n
        tid = threadid()
        for j in (i+1):n
            if A[i, j] > 0  # Edge exists
                # Use pre-computed neighbor lists
                neighbors_i = neighbors_list[i]
                neighbors_j = neighbors_list[j]

                common = intersect(neighbors_i, neighbors_j)
                for k in common
                    if k > j  # Maintain ordering i < j < k
                        push!(triangles[tid], (vertices[i], vertices[j], vertices[k]))
                    end
                end
            end
        end
    end
    
    # Flatten results and sort for deterministic output
    all_triangles = vcat(triangles...)
    sort!(all_triangles)
    return all_triangles
end

function detect_4_cycles(adj)
   
    vertices = sort(collect(keys(adj)))
    vertex_to_idx = Dict(v => i for (i, v) in enumerate(vertices))
    n = length(vertices)

    # Convert to sets for fast intersection
    adj_sets = Dict(v => Set(adj[v]) for v in vertices)

    # Early pruning: filter out vertices with degree < 2 (cannot participate in 4-cycles)
    vertices_deg2plus = filter(v -> length(adj_sets[v]) >= 2, vertices)
    n_filtered = length(vertices_deg2plus)
    
    # Use thread-local sets to store unique 4-cycles (represented in canonical form)
    cycles_per_thread = Vector{Set{Tuple{Any,Any,Any,Any}}}(undef, nthreads())
    for i in 1:nthreads()
        cycles_per_thread[i] = Set{Tuple{Any,Any,Any,Any}}()
    end

    # Find 4-cycles without materializing all wedges
    # For each pair of vertices (u, w), find their common neighbors
    vertex_pairs = [(vertices_deg2plus[i], vertices_deg2plus[j]) for i in 1:n_filtered for j in (i+1):n_filtered]

    @threads for (u, w) in vertex_pairs
        tid = threadid()
        # Find common neighbors of u and w (these form the centers of wedges)
        if haskey(adj_sets, u) && haskey(adj_sets, w)
            centers = collect(intersect(adj_sets[u], adj_sets[w]))

            if length(centers) >= 2
                # Each pair of centers forms a 4-cycle with endpoints u, w
                for i in 1:length(centers)
                    for j in (i+1):length(centers)
                        v1, v2 = centers[i], centers[j]
                        # Store 4-cycle in actual cyclic order: u -> v1 -> w -> v2 -> u
                        # To ensure uniqueness, normalize to start with smallest vertex
                        cycle_order = [u, v1, w, v2]
                        min_idx = argmin(cycle_order)
                        # Rotate to start with minimum vertex
                        normalized_cycle = [cycle_order[mod1(min_idx + i - 1, 4)] for i in 1:4]
                        # Normalize direction: ensure second element < last element to avoid counting both directions
                        if normalized_cycle[2] > normalized_cycle[4]
                            # Reverse the cycle (keeping first element as anchor)
                            normalized_cycle = [normalized_cycle[1], normalized_cycle[4], normalized_cycle[3], normalized_cycle[2]]
                        end
                        canonical_cycle = tuple(normalized_cycle...)

                        push!(cycles_per_thread[tid], canonical_cycle)
                    end
                end
            end
        end
    end

    # Merge thread-local sets
    cycles_set = Set{Tuple{Any,Any,Any,Any}}()
    for thread_set in cycles_per_thread
        union!(cycles_set, thread_set)
    end

    # Convert set to sorted array for deterministic output
    all_cycles = collect(cycles_set)
    sort!(all_cycles)
    return all_cycles
end
end