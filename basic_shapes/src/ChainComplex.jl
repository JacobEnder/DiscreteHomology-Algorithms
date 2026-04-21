#===================================================================#
# ChainComplex.jl
# Build chain groups C_3 -> C_2 -> C_1 -> C_0 and boundary matrices
# over Z/2, using a finite catalog of Q^3 quotient templates.
#
# C_0: vertices of G
# C_1: edges of G (no loops, normalized i<j)
# C_2: all triangles + all 4-cycles in G (normalized)
# C_3: boundaries pushed forward from template embeddings
#
# IMPORTANT: For H_2 = ker(∂_2)/im(∂_3), only im(∂_3) matters.
# We therefore build C_3 as a spanning set for im(∂_3) by
# deduplicating identical boundary columns.
#===================================================================#

using SparseArrays

"""Container holding the chain complex data used to compute H_2."""
struct ChainComplex
    n_vertices::Int
    edges::Vector{Tuple{Int,Int}}
    n_edges::Int
    triangles::Vector{NTuple{3,Int}}
    four_cycles::Vector{NTuple{4,Int}}
    n_2cells::Int
    c3_supports::Vector{Vector{Int}}   # each column: list of C2 indices with 1s
    n_3cells::Int
    d1::SparseMatrixCSC{Bool, Int}
    d2::SparseMatrixCSC{Bool, Int}
    d3::SparseMatrixCSC{Bool, Int}
    edge_to_idx::Dict{Tuple{Int,Int}, Int}
    cycle_to_idx::Dict{Any, Int}
end

"""Edges of a cycle as normalized (i,j) with i<j."""
function get_cycle_edges(cycle::NTuple{N,Int}) where {N}
    es = Tuple{Int,Int}[]
    for i in 1:N
        v1 = cycle[i]
        v2 = cycle[(i % N) + 1]
        push!(es, v1 < v2 ? (v1, v2) : (v2, v1))
    end
    return es
end

"""Build ∂_1: C_1 -> C_0 (edge to its two endpoints)."""
function build_d1(n_vertices::Int, edges::Vector{Tuple{Int,Int}})
    I, J, V = Int[], Int[], Bool[]
    for (col, (u, v)) in enumerate(edges)
        push!(I, u); push!(J, col); push!(V, true)
        push!(I, v); push!(J, col); push!(V, true)
    end
    return sparse(I, J, V, n_vertices, length(edges))
end

"""Build ∂_2: C_2 -> C_1 from triangles and 4-cycles."""
function build_d2(edge_to_idx::Dict{Tuple{Int,Int},Int},
                  triangles::Vector{NTuple{3,Int}},
                  four_cycles::Vector{NTuple{4,Int}})
    n_edges = length(edge_to_idx)
    n_tri = length(triangles)
    n_sq = length(four_cycles)
    n_2 = n_tri + n_sq

    I, J, V = Int[], Int[], Bool[]

    for (j, tri) in enumerate(triangles)
        for e in get_cycle_edges(tri)
            idx = get(edge_to_idx, e, 0)
            idx == 0 && continue
            push!(I, idx); push!(J, j); push!(V, true)
        end
    end
    for (k, sq) in enumerate(four_cycles)
        col = n_tri + k
        for e in get_cycle_edges(sq)
            idx = get(edge_to_idx, e, 0)
            idx == 0 && continue
            push!(I, idx); push!(J, col); push!(V, true)
        end
    end

    return n_2 == 0 ? sparse(Int[], Int[], Bool[], n_edges, 0) : sparse(I, J, V, n_edges, n_2)
end

"""Build C_3 supports by pushing forward template boundaries through embeddings."""
function build_c3_supports(G::SimpleGraph,
                           cycle_to_idx::Dict{Any,Int},
                           templates::Vector{Template};
                           verbose::Bool=false)
    supports = Vector{Vector{Int}}()
    seen = Set{String}()

    for t in templates
        emb_maps = find_injective_homs(t.graph, G)
        isempty(emb_maps) && continue
        for f in emb_maps
            counts = Dict{Int,Int}()

            for cyc in t.boundary
                mapped = if cyc isa NTuple{3,Int}
                    (f[cyc[1]], f[cyc[2]], f[cyc[3]])
                else
                    cyc4 = cyc::NTuple{4,Int}
                    (f[cyc4[1]], f[cyc4[2]], f[cyc4[3]], f[cyc4[4]])
                end
                norm = normalize_cycle(mapped)
                idx = get(cycle_to_idx, norm, 0)
                if idx == 0
                    verbose && println("    [warn] mapped face $norm not present in C_2(G)")
                    continue
                end
                counts[idx] = get(counts, idx, 0) + 1
            end

            supp = Int[]
            for (idx, c) in counts
                isodd(c) && push!(supp, idx)
            end
            isempty(supp) && continue
            sort!(supp)

            key = join(supp, ",")
            key in seen && continue
            push!(seen, key)
            push!(supports, supp)
        end
    end

    return supports
end

"""Build ∂_3: C_3 -> C_2 from supports."""
function build_d3(n_2cells::Int, supports::Vector{Vector{Int}})
    I, J, V = Int[], Int[], Bool[]
    for (col, supp) in enumerate(supports)
        for row in supp
            push!(I, row); push!(J, col); push!(V, true)
        end
    end
    n_3 = length(supports)
    return n_3 == 0 ? sparse(Int[], Int[], Bool[], n_2cells, 0) : sparse(I, J, V, n_2cells, n_3)
end

"""Verify that dA*dB = 0 over Z/2 by explicit integer multiplication."""
function verify_chain_zero(dA::SparseMatrixCSC{Bool,Int}, dB::SparseMatrixCSC{Bool,Int}; label::String="")
    (size(dA, 2) == size(dB, 1)) || error("matrix sizes incompatible for chain check")
    A = Int.(Matrix(dA))
    B = Int.(Matrix(dB))
    P = mod.(A * B, 2)
    any(P .!= 0) && error("Chain condition failed: $label has $(count(P .!= 0)) nonzero entries")
    return true
end

"""
    build_chain_complex(G; templates, verbose)

Build the chain complex using the given template catalog.
"""
function build_chain_complex(G::SimpleGraph; templates::Vector{Template}, verbose::Bool=true)
    n_vertices = length(G.vertices)

    # C1
    edges = non_loop_edges(G)
    n_edges = length(edges)
    edge_to_idx = Dict{Tuple{Int,Int},Int}()
    for (i, e) in enumerate(edges)
        edge_to_idx[e] = i
    end

    verbose && println("  C_0: $n_vertices vertices")
    verbose && println("  C_1: $n_edges edges")

    # C2
    adj = adjacency_no_loops(G)
    verbose && println("  Detecting 3-cycles and 4-cycles...")
    triangles = detect_3_cycles(adj, n_vertices)
    four_cycles = detect_4_cycles(adj, n_vertices)
    n_2cells = length(triangles) + length(four_cycles)
    verbose && println("  C_2: $n_2cells 2-cells ($(length(triangles)) triangles + $(length(four_cycles)) squares)")

    cycle_to_idx = build_cycle_to_idx(triangles, four_cycles)

    # d1, d2
    d1 = build_d1(n_vertices, edges)
    d2 = build_d2(edge_to_idx, triangles, four_cycles)

    # C3 supports and d3
    verbose && println("  Building C_3 from template embeddings (deduplicated by boundary column)...")
    supports = build_c3_supports(G, cycle_to_idx, templates; verbose=false)
    n_3cells = length(supports)
    verbose && println("  C_3: $n_3cells columns spanning im(∂_3)")
    d3 = build_d3(n_2cells, supports)

    # Chain checks
    if n_edges > 0 && n_2cells > 0
        verify_chain_zero(d1, d2; label="∂_1∘∂_2")
        verbose && println("  ✓ Verified: ∂_1 ∘ ∂_2 = 0 (mod 2)")
    end
    if n_2cells > 0 && n_3cells > 0
        verify_chain_zero(d2, d3; label="∂_2∘∂_3")
        verbose && println("  ✓ Verified: ∂_2 ∘ ∂_3 = 0 (mod 2)")
    end

    return ChainComplex(
        n_vertices,
        edges, n_edges,
        triangles, four_cycles, n_2cells,
        supports, n_3cells,
        d1, d2, d3,
        edge_to_idx, cycle_to_idx,
    )
end
