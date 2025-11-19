module Homology

using SparseArrays

export add_columns!, compute_h0!, compute_h1!, compute_homology!, rank_mod2

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


function rank_mod2!(matrix::SparseMatrixCSC{Bool, Int})
    n_rows, n_cols = size(matrix)

    # Convert to BitMatrix for efficient bitwise operations
    bit_matrix = BitMatrix(undef, n_rows, n_cols)

    # Sparse-to-bit conversion so that we can treat rows and cols as BitVectors
    for col in 1:n_cols
        bit_matrix[:, col] .= false
        cs, ce = matrix.colptr[col], matrix.colptr[col+1] - 1
        for idx in cs:ce
            if matrix.nzval[idx]
                bit_matrix[matrix.rowval[idx], col] = true
            end
        end
    end

    rank = 0
    current_row = 1

    # Gaussian elimination with bit operations
    for col in 1:n_cols
        # Find pivot using bit scanning
        pivot_row = 0
        for row in current_row:n_rows
            if bit_matrix[row, col]
                pivot_row = row
                break
            end
        end

        if pivot_row == 0
            continue
        end

        # Row swap using BitVector operations
        if pivot_row != current_row
            temp = bit_matrix[current_row, :]
            bit_matrix[current_row, :] = bit_matrix[pivot_row, :]
            bit_matrix[pivot_row, :] = temp
        end

        # Elimination using vectorized XOR (operates on 64-bit chunks)
        pivot_row_view = view(bit_matrix, current_row, :)
        for row in 1:n_rows
            if row != current_row && bit_matrix[row, col]
                view(bit_matrix, row, :) .⊻= pivot_row_view
            end
        end

        rank += 1
        current_row += 1

        if current_row > n_rows
            break
        end
    end

    return rank
end

function compute_h0!(matrix::SparseMatrixCSC{Bool, Int}, structure_info; rank::Bool=false)

    # dim(H0) = |G_V| - rank(M_V).
    vertex_count = length(structure_info.vertex_cols)
    vertex_matrix = matrix[structure_info.vertex_rows, structure_info.edge_cols]
    boundary_rank = rank_mod2!(copy(vertex_matrix))
    return vertex_count - boundary_rank
   
end


function compute_h1!(matrix::SparseMatrixCSC{Bool, Int}, structure_info; rank::Bool=false)

    # From the paper, dim(H1) = |G_E| - rank(M_E) - rank(M_C). We can look 
    edge_count = length(structure_info.edge_cols)

    edge_matrix = matrix[structure_info.vertex_rows, structure_info.edge_cols]
    edge_boundary_rank = rank_mod2!(copy(edge_matrix))

    triangle_matrix = matrix[structure_info.edge_rows, structure_info.short_cycle_cols]
    triangle_boundary_rank = rank_mod2!(copy(triangle_matrix))

    return edge_count - edge_boundary_rank - triangle_boundary_rank

end

"""
Modal homology computation, so that we can use persistence or direct rank computation.
"""
function compute_homology!(matrix::SparseMatrixCSC{Bool, Int}, structure_info; rank::Bool=false)

    h0 = compute_h0!(copy(matrix), structure_info; rank=true)
    h1 = compute_h1!(copy(matrix), structure_info; rank=true)
    return (h0, h1) 

end


end # module