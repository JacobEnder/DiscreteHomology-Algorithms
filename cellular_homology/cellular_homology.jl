include("boundary.jl")
include("homology.jl")
using .Homology
using .Boundary
using JSON
using SparseArrays

"""
Read and parse graph data from JSON file
"""
function read_graph_data(filename::String)
    filepath = joinpath("..", "data", filename)
    
    if !isfile(filepath)
        error("File not found: $filepath")
    end
    
    data = JSON.parsefile(filepath)
    return data
end

"""
Process graphs from JSON file
"""
function process_graphs_from_file(filename::String)
    graphs_data = read_graph_data(filename)
    
    results = []

    println("Cellular homology")
    
    for graph_data in graphs_data
        println("\n" * "="^50)
        name = graph_data["name"]
        println("Processing graph: $name")

        boundary_matrix, structure = construct_boundary_matrix(graph_data)
        
        # Compute homology.
        homology_time = @elapsed h1 = Homology.compute_h1!(boundary_matrix, structure)
       
        println("Homology computation time: $(homology_time) seconds")
        
        #println("H0 Dimension: $h0")
        println("H1 Dimension: $h1")
        
    end
    
    return results
end

function main()
    
    filename = "std_var_edge_high_verts_4.json"
    
    try
        results = process_graphs_from_file(filename)

        return results
        
    catch e
        println("Error: $e")
        return nothing
    end
end

main()