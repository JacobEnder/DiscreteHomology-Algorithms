#===================================================================#
# TemplateMinimize.jl
#
# Greedy reduction of a template catalog to a smaller universal
# generating family.
#
# A template t=(X,beta) is redundant relative to a set S if beta lies
# in the span of pushforwards of templates in S along injective
# homomorphisms into X. This is a finite check done on the small graph
# X itself.
#
# This file is optional; the main pipeline works with any template list.
#===================================================================#

"""C_2 basis (triangles + squares) and lookup for a small graph."""
function c2_basis(G::SimpleGraph)
    n = length(G.vertices)
    adj = adjacency_no_loops(G)
    tris = detect_3_cycles(adj, n)
    sqs = detect_4_cycles(adj, n)
    dict = build_cycle_to_idx(tris, sqs)
    return (tris, sqs, dict)
end

"""Support (C2 indices) of a boundary cycle list inside a specific graph basis."""
function boundary_support(boundary::Vector{Union{NTuple{3,Int},NTuple{4,Int}}},
                          cycle_to_idx::Dict{Any,Int})
    counts = Dict{Int,Int}()
    for cyc in boundary
        idx = get(cycle_to_idx, cyc, 0)
        idx == 0 && continue
        counts[idx] = get(counts, idx, 0) + 1
    end
    supp = [idx for (idx, c) in counts if isodd(c)]
    sort!(supp)
    return supp
end

"""Convert support vectors to a BitMatrix columns container."""
function supports_to_bitmatrix(n_rows::Int, supports::Vector{Vector{Int}})
    n_cols = length(supports)
    A = BitMatrix(falses(n_rows, n_cols))
    for (j, supp) in enumerate(supports)
        for i in supp
            1 <= i <= n_rows && (A[i, j] = true)
        end
    end
    return A
end

"""Rank over GF(2) for a BitMatrix (in-place row reduction on a copy)."""
function rank_gf2_bit(A::BitMatrix)
    n_rows, n_cols = size(A)
    r = 1
    B = copy(A)
    for c in 1:n_cols
        pivot = 0
        for i in r:n_rows
            if B[i, c]
                pivot = i
                break
            end
        end
        pivot == 0 && continue
        if pivot != r
            B[r, :], B[pivot, :] = B[pivot, :], B[r, :]
        end
        for i in 1:n_rows
            (i == r || !B[i, c]) && continue
            B[i, :] .⊻= B[r, :]
        end
        r += 1
        r > n_rows && break
    end
    return r - 1
end

"""Compute all pushed-forward boundary supports of template s inside target graph Xt."""
function pushed_supports(s::Template, Xt::SimpleGraph, cycle_to_idx_t::Dict{Any,Int})
    outs = Vector{Vector{Int}}()
    maps = find_injective_homs(s.graph, Xt)
    isempty(maps) && return outs

    for f in maps
        counts = Dict{Int,Int}()
        for cyc in s.boundary
            mapped = if cyc isa NTuple{3,Int}
                (f[cyc[1]], f[cyc[2]], f[cyc[3]])
            else
                cyc4 = cyc::NTuple{4,Int}
                (f[cyc4[1]], f[cyc4[2]], f[cyc4[3]], f[cyc4[4]])
            end
            norm = normalize_cycle(mapped)
            idx = get(cycle_to_idx_t, norm, 0)
            idx == 0 && continue
            counts[idx] = get(counts, idx, 0) + 1
        end
        supp = [idx for (idx, c) in counts if isodd(c)]
        isempty(supp) && continue
        sort!(supp)
        push!(outs, supp)
    end
    return outs
end

"""Return true iff t.boundary is in the span of pushforwards from selected into t.graph."""
function boundary_in_span(t::Template, selected::Vector{Template})
    isempty(selected) && return false
    _, _, cycle_to_idx_t = c2_basis(t.graph)
    n_rows = length(cycle_to_idx_t)

    # Collect supports of pushed boundaries
    cols = Vector{Vector{Int}}()
    seen = Set{String}()
    for s in selected
        for supp in pushed_supports(s, t.graph, cycle_to_idx_t)
            key = join(supp, ",")
            key in seen && continue
            push!(seen, key)
            push!(cols, supp)
        end
    end
    isempty(cols) && return false

    # Rank test: does adding t's boundary change rank?
    A = supports_to_bitmatrix(n_rows, cols)
    r0 = rank_gf2_bit(A)

    b_supp = boundary_support(t.boundary, cycle_to_idx_t)
    isempty(b_supp) && return true
    Ab = BitMatrix(hcat(A, supports_to_bitmatrix(n_rows, [b_supp])))
    r1 = rank_gf2_bit(Ab)
    return r1 == r0
end

"""Greedy minimization of a template catalog (sound universal redundancy elimination)."""
function minimize_template_catalog(templates::Vector{Template})
    # sort small-first to maximize chance embeddings exist
    ts = sort(templates; by = t -> (length(t.graph.vertices), length(t.boundary)))
    selected = Template[]
    for t in ts
        isempty(t.boundary) && continue
        boundary_in_span(t, selected) && continue
        push!(selected, t)
    end
    return selected
end
