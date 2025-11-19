module Boundary

using JSON
using SparseArrays
using Base.Threads

include("cycle_detection.jl")
using .CycleDetection

export extract_edges, construct_boundary_matrix, get_cycle_edges

"""
Extract edges as ordered pairs from adjacency list representation
Returns a vector of edges as tuples (u, v) where u < v, using integer vertex IDs
Also returns vertex mapping dictionaries for conversion between string and integer IDs
"""
function extract_edges(adj_list::Dict)

    # Create integer vertex mapping for faster comparisons
    vertices = sort(collect(keys(adj_list)))
    vertex_to_int = Dict(v => i for (i, v) in enumerate(vertices))
    int_to_vertex = Dict(i => v for (i, v) in enumerate(vertices))
    
    # BitSets are faster than Sets (we're working mod 2)
    n_vertices = length(vertices)
    edge_bits = BitSet()
    edges = Vector{Tuple{Int, Int}}()
    
    # Estimate number of edges (this is just for a little speed boost)
    estimated_edges = sum(length(neighbors) for neighbors in values(adj_list)) ÷ 2
    sizehint!(edges, estimated_edges)
    
    for (vertex_str, neighbors) in adj_list
        vertex_int = vertex_to_int[vertex_str]
        
        for neighbor_str in neighbors
            if haskey(vertex_to_int, neighbor_str)  
                neighbor_int = vertex_to_int[neighbor_str]
                
                # Ensure consistent ordering: smaller vertex first
                u, v = vertex_int < neighbor_int ? (vertex_int, neighbor_int) : (neighbor_int, vertex_int)
                
                edge_code = u * n_vertices + v
                
                # Only add if not seen before 
                if edge_code ∉ edge_bits
                    push!(edge_bits, edge_code)
                    push!(edges, (u, v))
                end
            end
        end
    end
    
    # Sort by first vertex, then second vertex
    sort!(edges)
    
    return edges, vertex_to_int, int_to_vertex
end

# Extract constituent edges from a cycle
function get_cycle_edges(cycle, vertex_to_int::Dict)

    # Handle 3-cycles
    if length(cycle) == 3

        v1, v2, v3 = cycle  # sorted lexicographically
        # Convert to integers and ensure ordering
        i1, i2, i3 = vertex_to_int[v1], vertex_to_int[v2], vertex_to_int[v3]
        edges = [(i1, i2), (i1, i3), (i2, i3)]

        # Normalize each edge (smaller vertex first)
        return [(min(u, v), max(u, v)) for (u, v) in edges]

    # Handle 4-cycles
    elseif length(cycle) == 4

        # 4-cycle: vertices are in actual cyclic order
        v1, v2, v3, v4 = cycle
        # Convert to integers to preserve ordering
        i1, i2, i3, i4 = vertex_to_int[v1], vertex_to_int[v2], vertex_to_int[v3], vertex_to_int[v4]
        edges = [(i1, i2), (i2, i3), (i3, i4), (i4, i1)]

        # Normalize each edge (smaller vertex first)
        return [(min(u, v), max(u, v)) for (u, v) in edges]
    end
end

# Given a graph G, build its boundary matrix
function construct_boundary_matrix(graph_data::Dict)

    # Pull the vertices and adjacency list
    vertices = sort(graph_data["vertices"])
    adj_list = graph_data["adjacency_list"]

    # Extract edges with integer mapping (already sorted)
    edges, vertex_to_int, int_to_vertex = extract_edges(adj_list)

    # Detect short cycles (3-cycles and 4-cycles)
    triangles = detect_3_cycles(adj_list)
    four_cycles = detect_4_cycles(adj_list)
    num_four_cycles = length(four_cycles)
    short_cycles = vcat(triangles, four_cycles)

    # Create mappings for indexing (edges are now integer tuples)
    edge_to_idx = Dict(e => i for (i, e) in enumerate(edges))

    # Dimensions
    n_vertices = length(vertices)
    n_edges = length(edges)
    n_short_cycles = length(short_cycles)

    n_rows = n_vertices + n_edges
    n_cols = n_vertices + n_edges + n_short_cycles

    # Pre-allocate single vectors with estimated size to reduce memory load
    estimated_nnz = 2 * n_edges + 3 * n_short_cycles + n_short_cycles  # conservative estimate
    I = Vector{Int}()
    J = Vector{Int}()
    V = Vector{Bool}()  # Changed to Bool for mod 2 arithmetic
    sizehint!(I, estimated_nnz)
    sizehint!(J, estimated_nnz)
    sizehint!(V, estimated_nnz)

    # Process edge columns (vertex columns are all zeros, so we skip them)
    edge_col_offset = n_vertices

    for edge_idx in 1:n_edges
        edge = edges[edge_idx]
        u, v = edge  # These are now integer vertex IDs

        col = edge_col_offset + edge_idx

        # Add 1's at the vertex rows corresponding to endpoints (u, v are already integer indices)
        push!(I, u)
        push!(J, col)
        push!(V, true)

        push!(I, v)
        push!(J, col)
        push!(V, true)
    end

    # Process short cycle columns - pre-normalize cycle edges.
    cycle_col_offset = n_vertices + n_edges

    for cycle_idx in 1:n_short_cycles
        cycle = short_cycles[cycle_idx]
        col = cycle_col_offset + cycle_idx

        # Get edges of the cycle (already normalized to integers)
        cycle_edges = get_cycle_edges(cycle, vertex_to_int)

        for normalized_edge in cycle_edges
            # Direct lookup - edge_to_idx should have all edges (no need to normalize again)
            edge_idx_lookup = get(edge_to_idx, normalized_edge, 0)
            if edge_idx_lookup > 0
                edge_row = n_vertices + edge_idx_lookup

                push!(I, edge_row)
                push!(J, col)
                push!(V, true)
            end
        end
    end

    # Create sparse matrix
    boundary_matrix = sparse(I, J, V, n_rows, n_cols)

    # Store structure information. This is needed to (quickly) keep track of which columns correspond to which structures,
    # which column offsets to use, and so on.
    structure_info = (
        n_vertices = n_vertices,
        n_edges = n_edges,
        n_short_cycles = n_short_cycles,
        vertex_cols = 1:n_vertices,
        edge_cols = (n_vertices+1):(n_vertices+n_edges),
        short_cycle_cols = (n_vertices+n_edges+1):(n_vertices+n_edges+n_short_cycles),
        vertex_rows = 1:n_vertices,
        edge_rows = (n_vertices+1):(n_vertices+n_edges),
        vertex_to_int = vertex_to_int,
        int_to_vertex = int_to_vertex,
        edges = edges,
        short_cycles = short_cycles
    )

    return boundary_matrix, structure_info
end

end # module