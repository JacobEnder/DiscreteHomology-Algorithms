#===================================================================#
# Templates.jl
#
# A "3-cell template" for the Q^3-based boundary operator is a pair
#   (X, beta)
# where
# - X is a quotient of Q^3 by identifying vertices with equal images,
# - beta is the induced 2-chain beta = q_*(∂Q^3) in C_2(X) (mod 2),
#   expressed as a set of (triangle/square) cycles with odd multiplicity.
#
# This file builds the finite catalog of all such templates (up to
# graph isomorphism carrying beta to beta).
#===================================================================#

using SparseArrays
using JSON3

"""Template = a quotient graph X plus a mod-2 boundary 2-chain."""
struct Template
    name::String
    graph::SimpleGraph
    boundary::Vector{Union{NTuple{3,Int},NTuple{4,Int}}}  # cycles with odd parity
end

"""Compute the quotient graph of Q^3 given a quotient map q::Vector{Int}."""
function quotient_graph_from_q3_map(q::Vector{Int})
    @assert length(q) == 8
    k = maximum(q)
    vertices = collect(1:k)
    edge_set = Set{Tuple{Int,Int}}()
    # Pull edges from Q3 (including loops from reflexivity)
    Q3 = build_q3()
    for (u, v) in Q3.edges
        qu, qv = q[u], q[v]
        push!(edge_set, (qu, qv))
    end
    # Ensure loops exist
    for v in vertices
        push!(edge_set, (v, v))
    end
    edges = sort!(collect(edge_set))
    adjacency = build_adjacency(vertices, edges)
    return SimpleGraph("quotient", vertices, edges, adjacency)
end

"""Push the cube boundary through a vertex map f (length 8), mod 2."""
function boundary_from_q3_map(f::Vector{Int})
    counts = Dict{Any, Int}()
    for face in Q3_FACES
        cyc = face_cycle_under_map(face, f)
        cyc === nothing && continue
        norm = normalize_cycle(cyc)
        counts[norm] = get(counts, norm, 0) + 1
    end
    boundary = Union{NTuple{3,Int},NTuple{4,Int}}[]
    for (cyc, c) in counts
        (c % 2 == 1) && push!(boundary, cyc)
    end
    sort!(boundary)
    return boundary
end

#-------------------------------------------------------------------#
# Set partitions of {1,...,n} (Bell number enumeration)
#-------------------------------------------------------------------#

"""Enumerate all set partitions of {1,...,n}. Each partition is a vector of blocks."""
function enumerate_set_partitions(n::Int)
    n < 0 && error("n must be >= 0")
    if n == 0
        return Vector{Vector{Vector{Int}}}([Vector{Vector{Int}}()])
    end
    parts_prev = enumerate_set_partitions(n - 1)
    out = Vector{Vector{Vector{Int}}}()
    for p in parts_prev
        # insert n into existing blocks
        for bi in 1:length(p)
            p2 = [copy(b) for b in p]
            push!(p2[bi], n)
            sort!(p2[bi])
            push!(out, p2)
        end
        # new block {n}
        p3 = [copy(b) for b in p]
        push!(p3, [n])
        push!(out, p3)
    end
    return out
end

"""Canonical quotient map q[1..n] from a set partition (blocks ordered by min element)."""
function quotient_map_from_partition(partition::Vector{Vector{Int}}, n::Int)
    blocks = sort(partition; by = b -> minimum(b))
    q = zeros(Int, n)
    for (lbl, block) in enumerate(blocks)
        for v in block
            q[v] = lbl
        end
    end
    return q
end

#-------------------------------------------------------------------#
# Template catalog construction
#-------------------------------------------------------------------#

"""Check whether two templates are isomorphic (graph iso carrying boundary to boundary)."""
function templates_isomorphic(t1::Template, t2::Template)
    (iso, perm) = are_isomorphic(t1.graph, t2.graph)
    iso || return false

    # perm maps vertices of t1.graph -> vertices of t2.graph
    b1 = Set(Union{NTuple{3,Int},NTuple{4,Int}}(permute_cycle(c, perm)) for c in t1.boundary)
    b2 = Set(t2.boundary)
    return b1 == b2
end

"""
    build_template_catalog(; drop_degenerate=true, drop_zero_boundary=true)

Enumerate all (X,beta) template types coming from Q^3 partitions.
Deduplicate up to template isomorphism.
"""
function build_template_catalog(; drop_degenerate::Bool=true, drop_zero_boundary::Bool=true)
    parts = enumerate_set_partitions(8)
    templates = Template[]
    for part in parts
        q = quotient_map_from_partition(part, 8)
        drop_degenerate && is_degenerate_3cube_map(q) && continue

        X = quotient_graph_from_q3_map(q)
        beta = boundary_from_q3_map(q)
        drop_zero_boundary && isempty(beta) && continue

        t = Template("", X, beta)
        # dedupe
        dup = false
        for u in templates
            if templates_isomorphic(t, u)
                dup = true
                break
            end
        end
        dup && continue

        push!(templates, t)
    end

    # give stable names
    for (i, t) in enumerate(templates)
        n = length(t.graph.vertices)
        m = length(non_loop_edges(t.graph))
        b = length(t.boundary)
        templates[i] = Template("T$(i)_n$(n)_m$(m)_b$(b)", t.graph, t.boundary)
    end
    return templates
end

#-------------------------------------------------------------------#
# JSON serialization helpers (optional but convenient)
#-------------------------------------------------------------------#

"""Convert a SimpleGraph to a JSON-serializable Dict."""
function graph_to_dict(G::SimpleGraph)
    return Dict(
        "name" => G.name,
        "n" => length(G.vertices),
        "edges" => [ [u, v] for (u, v) in G.edges ],
    )
end

"""Convert a Template to a JSON-serializable Dict."""
function template_to_dict(t::Template)
    return Dict(
        "name" => t.name,
        "graph" => graph_to_dict(t.graph),
        "boundary" => [ collect(c) for c in t.boundary ],
    )
end

"""Rebuild a SimpleGraph from a Dict."""
function dict_to_graph(d::Dict)
    n = Int(d[:n])
    vertices = collect(1:n)
    raw_edges = Tuple{Int,Int}[(Int(e[1]), Int(e[2])) for e in d[:edges]]
    adjacency = build_adjacency(vertices, raw_edges)
    return SimpleGraph(string(d[:name]), vertices, raw_edges, adjacency)
end

"""Rebuild a Template from a Dict."""
function dict_to_template(d::Dict)
    G = dict_to_graph(Dict(d[:graph]))
    b = Union{NTuple{3,Int},NTuple{4,Int}}[]
    for cyc in d[:boundary]
        if length(cyc) == 3
            push!(b, (Int(cyc[1]), Int(cyc[2]), Int(cyc[3])))
        elseif length(cyc) == 4
            push!(b, (Int(cyc[1]), Int(cyc[2]), Int(cyc[3]), Int(cyc[4])))
        else
            error("boundary cycle must have length 3 or 4")
        end
    end
    sort!(b)
    return Template(string(d[:name]), G, b)
end

"""Write a template catalog to a JSON file."""
function write_template_catalog(filename::String, templates::Vector{Template})
    data = [template_to_dict(t) for t in templates]
    open(filename, "w") do io
        JSON3.write(io, data)
    end
end

"""Read a template catalog from a JSON file."""
function read_template_catalog(filename::String)
    data = JSON3.read(read(filename, String))
    return [dict_to_template(Dict(d)) for d in data]
end