#===================================================================#
# ComputeH2.jl
# Compute H_2 = ker(∂_2)/im(∂_3) over GF(2).
#
# Pure Julia row/column elimination over GF(2).
#===================================================================#

using SparseArrays

"""Result of H_2 computation."""
struct H2Result
    h2::Int
    rank_d2::Int
    rank_d3::Int
    nullity_d2::Int
    n_1cells::Int
    n_2cells::Int
    n_3cells::Int
    generators::Vector{Vector{Int}}  # each generator as list of C2 indices
end

"""Rank of a Bool sparse matrix over GF(2) via row reduction."""
function rank_gf2(M::SparseMatrixCSC{Bool,Int})
    n_rows, n_cols = size(M)
    (n_rows == 0 || n_cols == 0) && return 0
    A = BitMatrix(Matrix(M))
    r = 1
    for c in 1:n_cols
        pivot = 0
        for i in r:n_rows
            if A[i, c]
                pivot = i
                break
            end
        end
        pivot == 0 && continue
        if pivot != r
            A[r, :], A[pivot, :] = A[pivot, :], A[r, :]
        end
        # eliminate
        for i in 1:n_rows
            (i == r || !A[i, c]) && continue
            A[i, :] .⊻= A[r, :]
        end
        r += 1
        r > n_rows && break
    end
    return r - 1
end

"""Compute a basis for ker(M) over GF(2). Returns BitVectors of length n_cols."""
function kernel_basis_gf2(M::SparseMatrixCSC{Bool,Int})
    n_rows, n_cols = size(M)
    n_cols == 0 && return BitVector[]
    if n_rows == 0
        basis = BitVector[]
        for i in 1:n_cols
            v = falses(n_cols)
            v[i] = true
            push!(basis, v)
        end
        return basis
    end

    A = BitMatrix(Matrix(M))
    pivot_col = zeros(Int, n_rows)
    pivot_row = zeros(Int, n_cols)

    r = 1
    for c in 1:n_cols
        pivot = 0
        for i in r:n_rows
            if A[i, c]
                pivot = i
                break
            end
        end
        pivot == 0 && continue
        if pivot != r
            A[r, :], A[pivot, :] = A[pivot, :], A[r, :]
        end
        pivot_col[r] = c
        pivot_row[c] = r
        for i in 1:n_rows
            (i == r || !A[i, c]) && continue
            A[i, :] .⊻= A[r, :]
        end
        r += 1
        r > n_rows && break
    end

    basis = BitVector[]
    for c in 1:n_cols
        pivot_row[c] != 0 && continue
        x = falses(n_cols)
        x[c] = true
        for i in 1:n_rows
            pc = pivot_col[i]
            pc == 0 && continue
            if A[i, c]
                x[pc] = true
            end
        end
        push!(basis, x)
    end
    return basis
end

"""Compute generators of ker(d2)/im(d3) using a column-space test."""
function h2_generators(d2::SparseMatrixCSC{Bool,Int}, d3::SparseMatrixCSC{Bool,Int})
    _, n_2 = size(d2)
    n_2 == 0 && return Vector{Int}[]

    K = kernel_basis_gf2(d2)
    isempty(K) && return Vector{Int}[]

    # Build augmented matrix [d3 | K] and perform left-to-right column reduction.
    n_3 = size(d3, 2)
    nK = length(K)
    A = BitMatrix(zeros(Bool, n_2, n_3 + nK))
    if n_3 > 0
        A[:, 1:n_3] = BitMatrix(Matrix(d3))
    end
    for (i, k) in enumerate(K)
        A[:, n_3 + i] = k
    end

    row_pivot = zeros(Int, n_2)
    is_pivot_col = falses(n_3 + nK)

    for col in 1:(n_3 + nK)
        pivot_row = 0
        for row in 1:n_2
            if A[row, col] && row_pivot[row] == 0
                pivot_row = row
                break
            end
        end
        pivot_row == 0 && continue
        is_pivot_col[col] = true
        row_pivot[pivot_row] = col
        for other in (col+1):(n_3 + nK)
            A[pivot_row, other] || continue
            A[:, other] .⊻= A[:, col]
        end
    end

    gens = Vector{Int}[]
    for i in 1:nK
        col = n_3 + i
        if is_pivot_col[col]
            supp = findall(K[i])
            !isempty(supp) && push!(gens, supp)
        end
    end
    return gens
end

"""Compute H_2 from a ChainComplex."""
function compute_h2(cc::ChainComplex)
    rank_d2 = cc.n_2cells > 0 ? rank_gf2(cc.d2) : 0
    rank_d3 = cc.n_3cells > 0 ? rank_gf2(cc.d3) : 0
    nullity_d2 = cc.n_2cells - rank_d2
    h2 = nullity_d2 - rank_d3
    h2 < 0 && (h2 = 0)
    gens = h2_generators(cc.d2, cc.d3)
    return H2Result(h2, rank_d2, rank_d3, nullity_d2,
                    cc.n_edges, cc.n_2cells, cc.n_3cells, gens)
end

"""Pretty-print a generator as sum of triangles/squares."""
function describe_generator(gen_indices::Vector{Int}, cc::ChainComplex)
    parts = String[]
    n_tri = length(cc.triangles)
    for idx in gen_indices
        if idx <= n_tri
            push!(parts, "△" * string(cc.triangles[idx]))
        else
            push!(parts, "□" * string(cc.four_cycles[idx - n_tri]))
        end
    end
    return join(parts, " + ")
end
