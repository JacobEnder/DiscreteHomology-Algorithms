# Necessary packages
using SparseArrays; # for sparse matrix and vecotr handling 
using LinearAlgebra; # used for rank computations 
using Base.Threads; # parallel computing

# Add JSON package for parsing graph files
using JSON

# Define Graph Structure
mutable struct graph
    # araray of vertex labels
    vertices::Array 

    # array of edges in the graph. Needs to be simple graph so it should be symmetric and reflexive 
    edges::Array 
end   


#===================================================================#
# Functions relevant to making inputting graphs easier
#===================================================================#

#= Function to make a set symmetric and reflexive, useful for edge sets
    Inputs: 
        1) the set A
    Outputs: 
        1) the symmetric and reflexiove version of A
    
    Notes: 
        
=#
function makeSymmetricReflexive(A)
    B=copy(A)
    for a in A
        if !( (a[2],a[1]) in B)
            push!(B,(a[2],a[1]) )
        end
        if !( (a[1],a[1]) in B)
            push!(B,(a[1],a[1]) )
        end
        if !( (a[2],a[2]) in B)
            push!(B,(a[2],a[2]) )
        end
    end
    return B
end

#= Function to quotient a graph by relating two vertices
    Inputs: 
        1) the graph G
        2) the first vertex
        3) the second vertex
    Outputs: 
        1) the graph with v1 = v2
    
    Notes: 
        
=#
 function quotient_vertices(G::graph,v1::Any,v2::Any)
    new_vert=[]
    new_edges=[]

    for v in G.vertices
        if !(v==v2)
            push!(new_vert,v)
        end
    end

    for i in G.edges
        if !(i[1]==v2 || i[2]==v2)
            push!(new_edges,i)
        end
        if i[1]==v2 && !(i[2]==v2) && !(i[2]==v1)
            push!(new_edges,(v1,i[2]))
        end
        if !(i[1]==v2) && i[2]==v2 && !(i[1]==v1)
            push!(new_edges,(i[1],v1))
        end
    end
    return graph(new_vert,new_edges)
end

#= Function to generate the suspension of a graph
    Inputs: 
        1) the graph G
        2) the length of suspension
    Outputs: 
        1) the suspension of the graph
    
    Notes: 
        
=#
 function suspension(G::graph,n::Int)
    n=n-1
    vertices=[]
    edges=[]
    for i = 1:n
        for j in G.vertices
            push!(vertices,(j,i))
        end
    end
    for i in vertices
        for j in vertices
            if (i[1]==j[1]) & (abs(i[2]-j[2])==1)
                push!(edges,(i,j))
                push!(edges,(j,i))
            end
            if (i[2]==j[2]) & ((i[1],j[1]) in G.edges)
                push!(edges,(i,j))
                push!(edges,(j,i))
            end
            if i[2]==n
                push!(edges,(i,'s'))
                push!(edges,('s',i))
            end
            if i[2]==1
                push!(edges,(i,'n'))
                push!(edges,('n',i))
            end
        end 
    end
    push!(vertices,'n')
    push!(vertices,'s') 

    push!(edges,('s','s')) 
    push!(edges,('n','n')) 

    susp=graph(vertices,unique!(edges))
    return susp
end

#= Function to generate the cartesian product of an array of arrays
    Inputs: 
        1) the array of arrays
    Outputs: 
        1) the cartesian product
    
    Notes: 
        
=#
 function cartesian_product(arrays) 
    # Base case: if arrays is empty, return an array with an empty tuple
    if isempty(arrays)
        return [()]
    end
    
    # Get the first array and the rest of the arrays
    first_array, rest_arrays = arrays[1], arrays[2:end]
    
    # Recursive call to get the Cartesian product of the rest of the arrays
    rest_product = cartesian_product(rest_arrays)
    
    # Combine each element of the first array with each tuple from the rest product
    result = [(x, rest...) for x in first_array for rest in rest_product]
    
    return result
end

#= Function to compute the box product of a list of graphs
    Inputs: 
        1) the list of graphs
    Outputs: 
        1) the graph that is a box product of all graphs in the list
    
    Notes: 
        
=#
function multBoxProd(graphs) 

    vert=cartesian_product([G.vertices for G in graphs])
    edge=[]
    
    for v in vert
        
        for w in vert
            tot = 0 # check if its an edge in the product
            
            for i = 1:length(v)
                
                if !(v[i] == w[i]) && ( (v[i],w[i]) in graphs[i].edges) #check if theye connected in g_i
                    tot += 1
                elseif !((v[i],w[i]) in graphs[i].edges)
                    tot += 2
                end
                
            end
            if tot <2
                push!(edge,(v,w))
            end
        end
    end

    return(graph(vert,edge))

end

#= Function to generate the complete graph on n vertices
    Inputs: 
        1) a number n
    Outputs: 
        1) the complete graph on n vertices
    
    Notes: 
        
=#
function completeGraph(n)
    v=[]
    e=[]
    for i = 1:n
        push!(v,i)
    end

    for w1 in v
        for w2 in v
            push!(e,(w1,w2))
        end
    end

    return graph(v,e)
end

#= Function to relabel the vertices of a graph as integers
    Inputs: 
        1) the graph G
    Outputs: 
        1) the relabeled graph
    
    Notes: 
        
=#
function relabel_vertices(G::graph)
    new_vert=[]
    
    transformDict = Dict{Any,Any}()
    for i in 1:length(G.vertices)
        push!(new_vert,i)
        transformDict[G.vertices[i]] = i
    end
    B = Channel{Any}(length(G.edges))
    @threads for e in G.edges
        put!(B,(transformDict[e[1]],transformDict[e[2]]))
    end
    close(B)
    new_edge=collect(B)

    return graph(new_vert,new_edge)
end

#===================================================
NAIVE ALGORITHM
===================================================#

#= Function to generate the cartesian product of an array of arrays
    Inputs: 
        1) the array of arrays
    Outputs: 
        1) the cartesian product
    
    Notes: 
        
=#
function cartesian_product(arrays) 
    # Base case: if arrays is empty, return an array with an empty tuple
    if isempty(arrays)
        return [()]
    end
    
    # Get the first array and the rest of the arrays
    first_array, rest_arrays = arrays[1], arrays[2:end]
    
    # Recursive call to get the Cartesian product of the rest of the arrays
    rest_product = cartesian_product(rest_arrays)
    
    # Combine each element of the first array with each tuple from the rest product
    result = [(x, rest...) for x in first_array for rest in rest_product]
    
    return result
end

#= Function to compute the box product of a list of graphs
    Inputs: 
        1) the list of graphs
    Outputs: 
        1) the graph that is a box product of all graphs in the list
    
    Notes: 
        
=#
function multBoxProd(graphs) 

    vert=cartesian_product([G.vertices for G in graphs])
    edge=[]
    
    for v in vert
        
        for w in vert
            tot = 0 # check if its an edge in the product
            
            for i = 1:length(v)
                
                if !(v[i] == w[i]) && ( (v[i],w[i]) in graphs[i].edges) #check if theye connected in g_i
                    tot += 1
                elseif !((v[i],w[i]) in graphs[i].edges)
                    tot += 2
                end
                
            end
            if tot <2
                push!(edge,(v,w))
            end
        end
    end

    return(graph(vert,edge))

end

# Generates n cube as a graph
function nCubeG(n)
    if n==0
        return graph([[]],[]) # for 1 cube faces
    end
    I_1 = graph([0,1],[(0,0),(1,1),(0,1),(1,0)])
    if n==1
        return graph([[0],[1]],[([0],[0]),([1],[1]),([0],[1]),([1],[0])])
    end
    lst = []
    for i=1:n
        push!(lst,I_1)
    end

    return multBoxProd(lst)
end

# compute faces 
function faceMap(cube,lcube,map,i,eps)
    if length(map)==2
        return [map[eps+1]]
    end
    face=[]
    for v in lcube.vertices
        v2 = new_tuple = (v[1:i-1]..., eps, v[i:end]...) # insert eps
        n = findfirst(==(v2),cube.vertices)
        push!(face,map[n])
    end
    return face
end

function isDegenerate(cube,lcube,f, n)
    for i = 1:n
        f1 = faceMap(cube,lcube,f,i,0)
        f2 = faceMap(cube,lcube,f,i,1)
        if f1==f2 
            return true
        end
    end
    return false
end
# return non degen cubes
function nonDegeneracies(cube,lcube,maps)
    n = Int(log2(length(cube.vertices)))
    nonDegen = []

    for f in maps
        if !(isDegenerate(cube,lcube,f,n))
            push!(nonDegen,f)
        end
    end
    return nonDegen
end

function sparse_col_concat(vectors)
    # Determine the dimensions of the resulting matrix
    n_rows = length(vectors[1])
    n_cols = length(vectors)
    
    # Initialize arrays to store the row indices, column indices, and values of the non-zero elements
    row_indices = Int[]
    col_indices = Int[]
    values = eltype(vectors[1])[]
    
    # Iterate through each vector and extract its non-zero elements
    for (j, vec) in enumerate(vectors)
        for i in 1:length(vec)
            if vec[i] != 0
                push!(row_indices, i)
                push!(col_indices, j)
                push!(values, vec[i])
            end
        end
    end
    
    # Construct the sparse matrix using the collected non-zero elements
    return sparse(row_indices, col_indices, values, n_rows, n_cols)
end

# Get the boundary sum as a coord vec
function getCoordVec(f,lcubes,hcube,lcube,n)
    vec = spzeros(length(lcubes))
    for i = 1:n
        # negative face
        f1 = faceMap(hcube,lcube,f,i,0)
        m = findfirst(==(f1), lcubes)
        if !(m==nothing)
            vec[m] += (-1)^i
        end

        # positive face
        f1 = faceMap(hcube,lcube,f,i,1)
        m = findfirst(==(f1), lcubes)
        if !(m==nothing)
            vec[m] += (-1)^(i+1)
        end
    end
    return vec
end

# convert to matrix
function delNmatrix(hcubes,lcubes,hcube,lcube)
    V = []
    n = Int(log2(length(hcube.vertices)))
    for f in hcubes 
        vec = getCoordVec(f,lcubes,hcube,lcube,n)
        push!(V,vec)
    end
    if length(V)==0
        return []
    end
    M=sparse_col_concat(V)
    return M
end


function rank_mod2(A)
    m, n = size(A)
    # Create a dense Boolean matrix B representing A mod 2.
    B = falses(m, n)
    for i in 1:m
        for j in 1:n
            # Since each entry is a rational number representing an integer,
            # we can convert it to an Int safely.
            # Set B[i, j] to true if A[i, j] is odd.
            if mod(Int(A[i, j]), 2) == 1
                B[i, j] = true
            end
        end
    end

    # Gaussian elimination over GF(2)
    rank = 0
    pivot_row = 1
    for j in 1:n
        # Find a pivot in column j starting at pivot_row.
        pivot = 0
        for i in pivot_row:m
            if B[i, j]
                pivot = i
                break
            end
        end
        if pivot != 0
            # Swap the pivot row into position if necessary.
            if pivot != pivot_row
                B[pivot_row, :], B[pivot, :] = copy(B[pivot, :]), copy(B[pivot_row, :])
            end
            # Eliminate entries below the pivot using XOR (addition in GF(2)).
            for i in pivot_row+1:m
                if B[i, j]
                    B[i, :] .= xor.(B[i, :], B[pivot_row, :])
                end
            end
            rank += 1
            pivot_row += 1
            if pivot_row > m
                break
            end
        end
    end
    return rank
end

# Function for determining if a given pair of n-1 cubes forms an n cube. Returns true/false 
function isPairNcube(f,g,nDict)
    map=true
    for i in 1:length(f)
        if !(f[i] in nDict[g[i]])
            map=false
            break
        end
    end
    return map
end

# generates maps from n cube to graph given the collection of n-1 cubes 
function graphMaps2(A,nDict) # A = array of the n-1 cubes, G = graph
    B = [] # output array of n cubes
    for f in A
        for g in A
            if isPairNcube(f,g,nDict) 
                push!(B,vcat(f,g))
            end
        end 
    end
    return(B)
end


# section 4 improvements 
#= function to generate the n-cube
    Inputs: 
        1) a natural number n

    Outputs: 
        1) n-cube as an array of two arrays, first being the vertices, second being the faces
    
    Notes: 
        - kind of messy but works 
=#
function nCube(n::Int)

    # handling for zero cube
    if n==0
        return [[(0)],[]]
    end

    map=[[[0],[1]],[[[[0]],[[1]]]]]

    if n==1
        return map
    end

    
    # iteratively define higher dimensional cubes 

    for i=2:n
        new_map = [[],[]]
        p2 = deepcopy(map)

        map4 = deepcopy(map)

        for x = 1:length(map[1])
            push!(map4[1][x],0) # origonal becomes 0 face, new becomes 1 face

            push!(p2[1][x],1)
        end

        map = deepcopy(map4)
        
        map3=deepcopy(map)

        for x = 1:length(map[2]) # loop thru pairs of faces
            for y = 1:length(map[2][x][1]) # loop thru vertices in faces
                push!(map3[2][x][1][y],0)
                push!(map3[2][x][2][y],0)

                push!(p2[2][x][1][y],1)
                push!(p2[2][x][2][y],1)
                
            end
        end
        
        map=deepcopy(map3)
        

        new_map[1]=vcat(map[1],p2[1])
        for x = 1:length(map[2])
            push!(new_map[2],[vcat(map[2][x][1],p2[2][x][1]),vcat(map[2][x][2],p2[2][x][2])])
        end

        push!(new_map[2],[deepcopy(map[1]),deepcopy(p2[1])])

        map = deepcopy(new_map)

    

    end

    return map
end

function faceList(n)
    cube = nCube(n) # initialize n cube
    facesList=[] # output array of faces
    for f in cube[2]
        cface1=[] # negative face
        cface2=[] # positive face
        for i in f[1]
            m = findfirst(==(i), cube[1])
            push!(cface1,m)
        end

        for i in f[2]
            m = findfirst(==(i), cube[1])
            push!(cface2,m)
        end
        push!(facesList,[cface1,cface2]) # add the pair to the list
    end
    return facesList
end

#= Function to compute the faces of a map
    Inputs: 
        1) the map as an array
        2) the facelist given from the faceList function
    Outputs: 
        1) the faces of the map, paired as [negativeFace, positiveFace]
    
    Notes: 
        
=#
function faces(map,facesList)
    face=[]
    for f in facesList
        m1=[map[i] for i in f[1]] # negative face
        m2=[map[j] for j in f[2]] # positive face
        push!(face,[m1,m2])
    end
    return face
end

#= Function to generate a dictionary of coordinates of maps
    Inputs: 
        1) array of equivalence classes

    Outputs: 
        1) the coordinate dictionary
    
    Notes: 
        
=#
function getCoordDict(lowerCubes)
    cdict = Dict{}()
    for i = 1:length(lowerCubes)
        cdict[lowerCubes[i]] = i
    end

    return cdict 
end

# generates maps from n cube to graph given the collection of n-1 cubes 
function graphMaps3(A,nDict) # A = array of the n-1 cubes, G = graph
    B = [] # output array of n cubes
    for f in A
        for g in A
            if isPairNcube(f[1],g[1],nDict) 
                deg=[]
                for i in f[2]
                    if i in g[2]
                        push!(deg,i)
                    end
                end
                if f[1]==g[1]
                    push!(deg,Int(log2(length(f[1]))+1))
                end
                push!(B,[vcat(f[1],g[1]),deg])
            end
        end 
    end
    return(B)
end

# Get the boundary sum as a coord vec
function getCoordVec(f,coordDict,facesList)
    vec = spzeros(length(keys(coordDict)))
    fc = faces(f,facesList)
    for i =1:length(fc) 
        # negative face
        f1 = fc[i][1]
        if f1 in keys(coordDict)
            m = coordDict[f1]
            vec[m] += (-1)^i
        end

        # positive face
        f1 = fc[i][2]
        if f1 in keys(coordDict)
            m = coordDict[f1]
            vec[m] += (-1)^(i+1)
        end
    end
    return vec
end

# return non degen cubes
function nonDegeneracies(maps)
    nonDegen = []

    for f in maps
        if length(f[2])==0
            push!(nonDegen,f[1])
        end
    end
    return nonDegen
end

# convert to matrix
function delNmatrix(hcubes,facesList,coordDict)
    V = []
    for f in hcubes 
        vec = getCoordVec(f,coordDict,facesList)
        push!(V,vec)
    end
    if length(V)==0
        return []
    end
    M=sparse_col_concat(V)
    return M
end

# final graph maps
function finalGraphMaps(Ln,coordDict,facesList, nDict)
    V=[]
    for f in Ln
        for g in Ln
            if isPairNcube(f[1],g[1],nDict)
                deg = false
                for i in f[2]
                    if i in g[2]
                        deg = true
                    end
                end

                if f[1]==g[1]
                    deg = true
                end

                if !(deg)
                    vec = getCoordVec(vcat(f[1],g[1]),coordDict,facesList)
                    push!(V,vec)
                end
            end
        end
    end
    M = sparse_col_concat(V)
    return M 
end

function cubicalHomology(G,n)

    #initialize nhood dict and faces dict
    faceDict=Dict{}()
    for i=0:n+1
        faceDict[i] = faceList(i)
    end

    nDict = Dict{}()
    for v in G.vertices
        nDict[v]=[]
    end

    for e in G.edges
        push!(nDict[e[1]],e[2])
    end

    mapsD=Dict{}()
    mapsD[0]=[[[v],[]] for v in G.vertices]

    for i = 1:n
        mapsD[i] = graphMaps3(mapsD[i-1],nDict)
        GC.gc()  # Force garbage collection after each map generation
    end
    if n==0
        lnondeg = []
        nondeg = [f[1] for f in mapsD[0]]
    elseif n==1
        lnondeg = [f[1] for f in mapsD[0]]
        nondeg = nonDegeneracies(mapsD[n])
    else
        lnondeg = nonDegeneracies(mapsD[n-1])
        nondeg = nonDegeneracies(mapsD[n])
    end
   
    lcoordDict = getCoordDict(lnondeg)
    coordDict = getCoordDict(nondeg)

    if n==0
        dimKer = length(nondeg)
    else
        M1 = delNmatrix(nondeg,faceDict[n],lcoordDict)
        dimKer=size(M1,2)-rank_mod2(Matrix(M1))
        M1 = nothing  # Clear reference to M1
        GC.gc()  # Force garbage collection after computing kernel dimension
    end

    M2 = finalGraphMaps(mapsD[n], coordDict, faceDict[n+1],nDict)

    dimIM = rank_mod2(Matrix(M2))

    M2 = nothing  # Clear reference to M2
    GC.gc()  # Force garbage collection after computing image dimension

    return dimKer-dimIM
end

"""
Read and parse graph data from JSON file
"""
function read_graph_data(filename::String)
    # Read from ../data/ directory
    filepath = joinpath("..", "data", filename)
    
    if !isfile(filepath)
        error("File not found: $filepath")
    end
    
    data = JSON.parsefile(filepath)
    return data
end

"""
Convert JSON graph data to graph struct
"""
function json_to_graph(graph_data)
    # Extract vertices - convert from strings to integers
    vertices = [parse(Int, v) for v in graph_data["vertices"]]
    
    # Extract edges from adjacency_list (which is a dictionary)
    adjacency_list = graph_data["adjacency_list"]
    edges = []
    
    # Convert adjacency list to edge tuples
    for (vertex_str, neighbors) in adjacency_list
        vertex = parse(Int, vertex_str)
        for neighbor_str in neighbors
            neighbor = parse(Int, neighbor_str)
            push!(edges, (vertex, neighbor))
        end
    end
    
    # Make edges symmetric and reflexive
    edges = makeSymmetricReflexive(edges)
    
    return graph(vertices, edges)
end

"""
Process graphs from JSON file
"""
function process_graphs_from_file(filename::String)
    graphs_data = read_graph_data(filename)

    results = []

    println("Cubical Homology (Discrete Morse Theory)")

    for graph_data in graphs_data
        println("\n" * "="^50)
        name = graph_data["name"]
        println("Processing graph: $name")

        # Convert JSON to graph struct
        G = json_to_graph(graph_data)

        # Compute homology with timing
        homology_time = @elapsed h1 = cubicalHomology(G, 1)

        println("Homology computation time: $(homology_time) seconds")
        println("H1 Dimension: $h1")

        push!(results, Dict("name" => name, "h1" => h1, "time" => homology_time))

        # Clear graph data and force garbage collection after each graph
        G = nothing
        graph_data = nothing
        GC.gc()
    end

    return results
end

function main()

    filename = "std_var_edge_high_verts_4.json"

    try
        # Force garbage collection at start
        GC.gc()

        results = process_graphs_from_file(filename)

        # Force final garbage collection
        GC.gc()

        return results

    catch e
        println("Error: $e")
        return nothing
    end
end

main()