module Homology

using SparseArrays
using Modulo2

export lowest_nonzero_entry, add_columns!, compute_h0!, compute_h1!, compute_homology!, rank_mod2

function lowest_nonzero_entry(matrix::SparseMatrixCSC{Bool, Int}, col::Int)
    if col < 1 || col > size(matrix, 2)
        return -1
    end

    # Get column range in sparse format
    cs = matrix.colptr[col]
    ce = matrix.colptr[col+1] - 1

    if cs > ce  # Empty column
        return -1
    end

    # Find the maximum row index with nonzero entry
    max_row = -1
    for idx in cs:ce
        if matrix.nzval[idx]  # Bool values are already mod 2
            row = matrix.rowval[idx]
            if row > max_row
                max_row = row
            end
        end
    end

    return max_row
end

"""
Add column j to column i in mod 2 arithmetic (i := i + j mod 2).
Modifies the matrix in place using XOR.
"""
function add_columns!(matrix::SparseMatrixCSC{Bool, Int}, i::Int, j::Int)
    if i == j || i < 1 || j < 1 || i > size(matrix, 2) || j > size(matrix, 2)
        return
    end

    # For Bool matrices, addition mod 2 is simply XOR
    matrix[:, i] = matrix[:, i] .⊻ matrix[:, j]
end


"""
Compute rank modulo 2 using the Modulo2 package.
Converts the sparse Bool matrix to ZZ2Matrix.
"""
function rank_mod2!(matrix::SparseMatrixCSC{Bool, Int})
    n_rows, n_cols = size(matrix)
    
    # Convert sparse Bool matrix to dense Bool array
    dense_matrix = Matrix{Bool}(matrix)
    
    # Convert to ZZ2Matrix
    zz2_matrix = ZZ2Matrix(dense_matrix)
    
    # Use Modulo2's efficient rank computation
    return rank(zz2_matrix)
end

function compute_h0!(matrix::SparseMatrixCSC{Bool, Int}, structure_info)

    # dim(H0) = |G_V| - rank(M_V).
    
    vertex_count = length(structure_info.vertex_cols)
    vertex_matrix = matrix[structure_info.vertex_rows, structure_info.edge_cols]
    boundary_rank = rank_mod2!(copy(vertex_matrix))
    return vertex_count - boundary_rank
end


function compute_h1!(matrix::SparseMatrixCSC{Bool, Int}, structure_info)

    # From the paper, dim(H1) = |G_E| - rank(M_E) - rank(M_C). We can look 
    
    edge_count = length(structure_info.edge_cols)

    edge_matrix = matrix[structure_info.vertex_rows, structure_info.edge_cols]
    edge_boundary_rank = rank_mod2!(copy(edge_matrix))

    triangle_matrix = matrix[structure_info.edge_rows, structure_info.short_cycle_cols]
    triangle_boundary_rank = rank_mod2!(copy(triangle_matrix))

    return edge_count - edge_boundary_rank - triangle_boundary_rank
    
end

function compute_homology!(matrix::SparseMatrixCSC{Bool, Int}, structure_info)
        h0 = compute_h0!(copy(matrix), structure_info; rank=true)
        h1 = compute_h1!(copy(matrix), structure_info; rank=true)
        return (h0, h1)
end

end # module