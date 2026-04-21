#===================================================================#
# GraphTypes.jl
# Minimal graph type + JSON loader utilities used throughout the v2
# discrete H_2 pipeline.
#
# Conventions
# - Vertices are 1-indexed Ints.
# - Graphs are treated as reflexive (self-loops present at every
#   vertex) and symmetric.
# - When we talk about edges in C_1 we drop self-loops and normalize
#   to (i,j) with i<j.
#===================================================================#

using JSON3

"""A small adjacency-list graph with a name."""
mutable struct SimpleGraph
    name::String
    vertices::Vector{Int}
    edges::Vector{Tuple{Int,Int}}          # includes (v,v) loops
    adjacency::Dict{Int, Vector{Int}}      # includes loops
end

"""Symmetric reflexive closure of a directed edge list."""
function make_symmetric_reflexive(edges::Vector{Tuple{Int,Int}}, vertices::Vector{Int})
    edge_set = Set{Tuple{Int,Int}}()
    for (u, v) in edges
        push!(edge_set, (u, v))
        push!(edge_set, (v, u))
    end
    for v in vertices
        push!(edge_set, (v, v))
    end
    return sort!(collect(edge_set))
end

"""Build an adjacency dict from (possibly looped) edges."""
function build_adjacency(vertices::Vector{Int}, edges::Vector{Tuple{Int,Int}})
    adj = Dict{Int, Vector{Int}}()
    for v in vertices
        adj[v] = Int[]
    end
    for (u, v) in edges
        push!(adj[u], v)
    end
    for v in vertices
        adj[v] = sort!(unique!(adj[v]))
    end
    return adj
end

"""Return undirected non-loop edges normalized as (i,j) with i<j."""
function non_loop_edges(G::SimpleGraph)
    edge_set = Set{Tuple{Int,Int}}()
    for (u, v) in G.edges
        u == v && continue
        push!(edge_set, u < v ? (u, v) : (v, u))
    end
    return sort!(collect(edge_set))
end

"""Adjacency without self-loops (used by cycle enumeration)."""
function adjacency_no_loops(G::SimpleGraph)
    adj = Dict{Int, Vector{Int}}()
    for v in G.vertices
        adj[v] = [u for u in G.adjacency[v] if u != v]
        sort!(adj[v])
    end
    return adj
end

"""Parse graphs from the json format used by the original pipeline."""
function parse_graph_json(filename::String)
    content = read(filename, String)
    data = JSON3.read(content)

    graphs = SimpleGraph[]
    for entry in data
        name = string(entry["name"])
        orig_verts = [string(v) for v in entry["vertices"]]
        n = length(orig_verts)

        # relabeling map: original string label -> integer 1..n
        label_map = Dict{String, Int}()
        for (i, v) in enumerate(orig_verts)
            label_map[v] = i
        end

        raw_edges = Tuple{Int,Int}[]
        adj_list = entry["adjacency_list"]
        for (v_str, neighbors) in pairs(adj_list)
            v = label_map[string(v_str)]
            for nb in neighbors
                w = label_map[string(nb)]
                push!(raw_edges, (v, w))
            end
        end

        vertices = collect(1:n)
        edges = make_symmetric_reflexive(raw_edges, vertices)
        adjacency = build_adjacency(vertices, edges)
        push!(graphs, SimpleGraph(name, vertices, edges, adjacency))
    end
    return graphs
end
