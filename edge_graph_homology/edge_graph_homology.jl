# Necessary packages -- run "using Pkg Pkg.add("name")" to add packages
using SparseArrays; # for sparse matrix and vecotr handling 
using LinearAlgebra; # used for rank computations 
using Base.Threads; # parallel computing
using Modulo2; # computing over Z2
using Mods; #computing over finite fields
using Primes; # determining which field to compute over
using BenchmarkTools; # for benchmarking 
using JSON3; # for reading JSON files

# structure for graphs 

# Define Graph Structure
mutable struct graph
    # araray of vertex labels
    vertices::Array 

    # array of edges in the graph. Needs to be simple graph so it should be symmetric and reflexive 
    edges::Array 
end  

# function to compute the rank of a sparse matrix mod p using REF
function ref_rank!(A::SparseMatrixCSC{Mod{p}}, inplace = false) where p

    if !inplace
        A = copy(A)
    end
    
    m, n = size(A)

    # quicker to have less rows
    if m > n 
        A = sparse(Matrix(transpose(A)))
        m, n = size(A)
    end

    lead = 1  # Track the lead position

    rank = 0 # rank is number of pivots found

    for r in 1:m
        if lead > n
            return rank
        end

        
        # Find the pivot row
        i = r
        while lead <= n && A[i, lead] == 0
            
            i += 1
            if i > m # column zero below r
                i = r
                lead += 1
                if lead > n # All columns exhausted
                    return rank 
                end
            end
        end

        # pivot was found, increase rank 

        rank+=1

        # Swap rows if necessary

        if i != r
            A[[r,i], :] .= A[[i,r], :]
        end


        # Normalize the pivot row
        pivot = A[r, lead]
        
        
        A[r, :] = A[r, :] ./ pivot  
        
        
        # Eliminate entries below the pivot
        for i = r+1:m
            if A[i, lead] != 0
                A[i, :] -= A[i, lead] * A[r, :]
            end
        end

        lead += 1
    end

    return rank
end

#===================================================================#
# Functions for preprocessing graphs
#===================================================================#

#= Function that removes a vertex from a graph and nhood dict
    Inputs: 
        1) the vertex to remove
        2) the graph G 
        3) the nhood dict of the graph
    Outputs: 
        1) the new graph
        2) the new nhood dict
    
    Notes: 
        
=#
function remove_vertex(v,G,nhoodDict)

    # remove from graph
    filter!(x -> x != v, G.vertices)
    filter!(t -> v ∉ t, G.edges)    

    pop!(nhoodDict,v)

    for key in keys(nhoodDict)
        filter!(x -> x != v, nhoodDict[key])
    end

    return G, nhoodDict
end

#= Function to detect and removes a vertex of degree n,  if none found returns the graph along with a false flag
    Inputs: 
        1) the graph G
        2) the nhood dict of the graph
    Outputs: 
        1) the graph with one vertex removed (if found)
        2) the nhood dict with one vertex removed (if found)
        3) true/false value if vertex was removed
    
    Notes: 
        
=#
function remove_vert_deg_n(G,nhoodDict)

    for key in keys(nhoodDict) 
        for v in filter(x -> x != key, nhoodDict[key]) # dont consider loops

            # check if v is connected to everything the key is
            conn = true
            for w in filter(x -> x != key, nhoodDict[key]) # still dont consider loops bc already know key ~ v
                if !(w in nhoodDict[v])
                    conn = false
                    break
                end
            end

            if conn == true
                G, nhooodDict = remove_vertex(key,G,nhoodDict)
                return G, nhoodDict, true
            end
        end
    end

    return G, nhoodDict, false # no case found
end

#= Function to preproccess a graph, making homology computations quicker
    Inputs: 
        1) the graph G
    Outputs: 
        1) the optimal graph to compute on given results about preprocessing
    
    Notes: 
        
=#
function preprocess_graph(G::graph)
    G2 = graph(deepcopy(G.vertices), deepcopy(G.edges))
    nhoodDict = get_nhood_dict(G2)
    flg = true

    while flg
        G2, nhoodDict, flg = remove_vert_deg_n(G2,nhoodDict)
    end

    return G2
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
    return graph(new_vert,unique!(new_edges))
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

# random graph
function random_graph(n_vert,av_edge)
    V = 1:n_vert
    E=[]
    for x in V 
        for i=1:av_edge
            y = rand(V)
            push!(E,(x,y))
        end
        push!(E,(x,x))
    end
    E = unique!(E)
    E = makeSymmetricReflexive(E)

    return graph(V,E)
end

#=======================================================================================
Functions for generating the equivalence classes according to the hyperoctahedrial group
========================================================================================#

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
 
#= function to generate the reversals as lists of -1 and 1
    Inputs: 
        1) a natural number n, representing the dimension you want to generate (length of list)

    Outputs: 
        1) An array of all reversals for that dimension
    
    Notes: 
        - 1 does nothing to a coordinate, -1 changes form 1 to 0 and 0 to 1 
=#
 function generate_reversals(n::Int)
    if n == 0
        return [[]] # empty array for 0th dimension
    else
        shorter_lists = generate_reversals(n - 1)
        return vcat([[1; l] for l in shorter_lists], [[-1; l] for l in shorter_lists])
    end
end

#= function to generate all permutations of a list (used for symmetries)
    Inputs: 
        1) a list 

    Outputs: 
        1) a list of all permutations 
    
    Notes: 
        
=#
 function permutations(lst)
    if length(lst) == 0
        return [[]]
    else
        result = []
        for i in 1:length(lst)
            # Take the current element
            current = lst[i]
            # Remaining list without the current element
            remaining = [lst[j] for j in 1:length(lst) if j != i]
            # Recursively find permutations of the remaining list
            for perm in permutations(remaining)
                push!(result, [current; perm])
            end
        end
        return result
    end
end

#= function to determine the sign of a permutation
    Inputs: 
        1) a permutation

    Outputs: 
        1) its sign 
    
    Notes: 
        
=#
 function sign_of_permutation(perm)
    n = length(perm)
    inversions = 0
    for i in 1:(n-1)
        for j in (i+1):n
            if perm[i] > perm[j]
                inversions += 1
            end
        end
    end
    return (-1)^inversions
end

#= function to generate all permutations of a list (used for symmetries)
    Inputs: 
        1) an integer n

    Outputs: 
        1) a list of all elements of the n-th hyperoctahedrial group, written as [symmetry, reversal, sign] 
    
    Notes: 
        
=#
 function hyperOctahedrial(n::Int)
    grp = []

    list = 1:n 
    # Generate all permutations
    symmetric_Gp = permutations(list) # lists represent image of ordered list
    reversals = generate_reversals(n)
    for g in symmetric_Gp
        for r in reversals

            sgn = prod(r)*sign_of_permutation(g)

            push!(grp,[g,r,sgn])
        end
    end

    return grp
end

#= function to permute the coordinates of a vertex in the n-cube according to a group element in the hyperOctahedrial group
    Inputs: 
        1) vertex in coordinate representation
        2) the group element

    Outputs: 
        1) the permuted vertex
    
    Notes: 
        
=#
 function permuteCoords(vertex, gpElet)
    
    pm1=[]
    for n = 1:length(vertex)
        push!(pm1,vertex[gpElet[1][n]])
    end

    for n = 1:length(vertex)
        pm1[n] = Int(mod(pm1[n] + (gpElet[2][n]+3)/2 , 2))
    end
    return pm1
end

#= Function that calculates the image of an n-cube under a group element
    Inputs: 
        1) then n-cube
        2) the group element

    Outputs: 
        1) the image
    
    Notes: 
        
=#
function calculateImageCube(cube,gpElet)
    cube2=deepcopy(cube) # dont change origonal cube
    image = []

    for n=1:length(cube2[1])
        push!(image, permuteCoords(cube2[1][n],gpElet))
    end

    return image

end

#= Function to take a cube and its image under a transformation, and give an array representing the map from the cube to its image
    Inputs: 
        1) the n-cube
        2) the image

    Outputs: 
        1) an array representing the map
    
    Notes: transformation[i] = vertex in image that cube.vertices[i] gets sent to
        
=#
 function transformationCoords(cube,image)
    transformation = []
    for n = 1:length(image)
        m = findfirst(==(image[n]), cube[1]) # coord of image in origonal cube
        push!(transformation,m)
    end
    return transformation
end

#= Function to generate the list of the image maps, how the faces get swapped, and the sign of the permutation for the entire n-th hyperoctahedrial group
    Inputs: 
        1) a natural number n

    Outputs: 
        1) list of maps as an array of arrays
    
    Notes: 
        - each inner array is an array of three arrays, first is the map, second is how the faces are mapped, and third is the sign of the permutation
        
=#
function generateEQClist(n)
    EQClist = []

    cube = nCube(n)
    grp = hyperOctahedrial(n)

    for g in grp # if parralelized, need first thing in list to be identity
        transformation = transformationCoords(cube,calculateImageCube(cube, g))
        push!(EQClist,[transformation,g[3]])
    end
    
    return EQClist
end

#= Function that generates the equivalence class of a map
    Inputs: 
        1) the map as an array
        2) the list of maps from the generateEQClsit function

    Outputs: 
        1) the equivalence class of the map
    
    Notes: 
        
=#
function generateEQC(map,EQClist)
    EQC = [] # output array of maps
    for el in EQClist
        cmap = [map[im] for im in el[1]] # new map under transformation
        push!(EQC,[cmap,el[2]])
    end        
    return EQC
end

#= Function that creates a list of the faces of a cube in terms of position in the origonal map
    Inputs: 
        1) a natural number n

    Outputs: 
        1) a list of the faces 
    
    Notes: 
        - faces are paied [A,B], where A is negative face and B is positive
        
=#
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

#===================================================================#
# Functions to compute the homology of a Graph 
#===================================================================#


#= Function to generate a dictionary of neighborhoods of a graph
    Inputs: 
        1) the graph G

    Outputs: 
        1) the neighborhood dictionary of the graph
    
    Notes: 
        - dictionary of the form a => A, where A is the set of vertices a is connected to
=#
function get_nhood_dict(G)
    # initialize nhood dict for checking maps
    nhoodDict = Dict{Any,Array{Any}}()
    for v in G.vertices
        nhoodDict[v] = []
    end

    for e in G.edges
        push!(nhoodDict[e[1]], e[2])
    end
    return nhoodDict
end

#= Function to determine if a pair of n-1 cubes forms an n cube by making f the n-th negative face and g the n-th positive face
    Inputs: 
        1) the first n-1 cube
        2) the second n-1 cube
        3) the neighborhood dictionary of the graph

    Outputs: 
        1) true/false value of if they form an n cube
    
    Notes: 
        
=#
 function isPairNcube(f,g,nhoodDict)
    map=true
    
    for i in 1:length(f[1]) # since maps are pairs [A,B,C], where A is the map, B is the degenerate coordinates, and C is the sign relative to the EQC rep
        if !(g[1][i] in nhoodDict[f[1][i]]) #check that pairs are connected
            map=false
            break
        end
    end
    return map
end

#= Function that generates n-cubes of a graph
    Inputs: 
        1) an array of the n-1 cubes
        2) the graph G
        3) the equivalence class list for the n-th dimension
        4) the neighborhood dictionary of the graph

    Outputs: 
        1) the n cubes
    
    Notes: 
        
=#
function graphMaps(A,G,EQClist,nhoodDict) # A = array of the n-1 cubes, G = graph
    
    # check if generating zero cubes
    if length(A)==0
        C=[]
        for i in G.vertices
            push!(C,[[[i],1]]) # last coord is automorphism sign wrt first cube in class
        end
        return C
    end

    B = Channel{Vector{Any}}(length(A)*length(A)*length(A[1])*length(A[1])) # output channel of maps
    
    @threads for i = 1:length(A) 
        for j = i:length(A)
            f=A[i]
            g=A[j]
            
            for t=1:length(f)
                if isPairNcube(f[1],g[t],nhoodDict) 
                    h = vcat(f[1][1],g[t][1])

                    eq = generateEQC(h,EQClist)
                    
                    put!(B,eq)
                    
                end
            end
        end 
    end
    close(B) # close the channel
    cubes = collect(B) # collect the cubes
    return(cubes)
end


#= Function that checks if a map is semi degenerate, i.e. if x~-x
    Inputs: 
        1) the map
        2) the n-th equivalence class list

    Outputs: 
        1) a true/false value of if its semi degenerate
    
    Notes: 
        
=#
function checkSemiDegen(map,EQClist)

    for g in EQClist
        if g[2]==-1 # check all maps with sign -1 in EQC
            if map == [map[im] for im in g[1]] # check that they're equal
                return true
            end
        end
    end
    return false
    
end

#= Function to remove semi degenerate maps
    Inputs: 
        1) an array of maps (here they are equivalence classes)

    Outputs: 
        1) the non semi degenerate maps
    
    Notes: 
        
=#
function semiNonDegen(maps)
    
    B=Channel{Vector{Any}}(length(maps))
    @threads for i in maps
        semiDegen=false
        for j in i
            if i[1][2]*j[2] == -1
                if i[1][1]==j[1] # check if semi degenerate
                    semiDegen=true
                    break
                end
            end
        end
        if !(semiDegen)
            put!(B,i) # if not semi degenerate, add to list
        end
    end
    close(B)
    semiNonDegen = collect(B)
    return semiNonDegen
end

#= Function to check if a map f appears in the equivalence class g
    Inputs: 
        1) the map f
        2) the equivalence class g

    Outputs: 
        1) a true/false value
    
    Notes: 
        
=#
function is_related(f,g) # f = map, g = equivalence class

    if !(sort(f[1])==sort(g[1][1])) # check if they have same vertices
        return false
    end
    
    for h in g
        if h[1]==f[1] # check if equal 
            return true
        end
    end

    return false
end

#= Function to remove duplicate equivalence classes
    Inputs: 
        1) the array of equivalence calsses of maps

    Outputs: 
        1) the array of unique equivalence classes
    
    Notes: 
        
=#
function remove_duplicates(maps) 
    B=Channel{Vector{Any}}(length(maps))
    @threads for i = 1:length(maps)
        k=maps[i]
        f = k[1]
        not_dupe = true
        for j = i+1:length(maps)
            g = maps[j]
            if is_related(f,g)
                not_dupe = false
                break
            end
        end

        if not_dupe
            put!(B,k)
        end
    end

    close(B)
    C = collect(B)
    return C
end

#= Function to generate a dictionary of coordinates of maps
    Inputs: 
        1) array of equivalence classes

    Outputs: 
        1) the coordinate dictionary
    
    Notes: 
        
=#
function coordDict(lowerCubes)
    cdict = Dict{}()
    for i = 1:length(lowerCubes)
        for j = 1:length(lowerCubes[i])
            cdict[lowerCubes[i][j][1]] = [i,lowerCubes[i][j][2]]
        end
    end

    return cdict 
end

#= Function converts an array of sparse vectors to a matrix where vecotrs are column vecs
    Inputs: 
        1) the array of vectors

    Outputs: 
        1) the matrix
    
    Notes: 
        
=#
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
        inds, vals = findnz(vec)
        cols = fill(Int64(j),length(inds))

        row_indices = vcat(row_indices,inds) 
        col_indices = vcat(col_indices, cols) 
        values = vcat(values, vals)

    end
    
    
    return sparse(row_indices, col_indices, values, n_rows, n_cols)
    
end

#= Function computes the boundary of a map and stores it as a sparse vector
    Inputs: 
        1) the map as an array
        2) the lower non degenerate cubes
        3) the n-dimensional face list

    Outputs: 
        1) the boundary as a coordinate vector, stored sparsly
    
    Notes: 
        
=#
function boundarySum(map,cdict,faceList)
    
    c_rows = []
    c_vals = Rational{Int}[]
    
    fc = faces(map[1],faceList) # generate faces

    for i = 1:length(fc) # loop through faces
        # negative face
        key = fc[i][1]
        if key in keys(cdict)
            indSgn=cdict[key]
            ind = indSgn[1]
            sgn = indSgn[2]

            pos = findfirst(x->x==ind,c_rows)
            if pos == nothing 
                push!(c_rows,ind)
                push!(c_vals, Rational{Int}((-1)^(i)*sgn))
            else
                c_vals[pos]+=Rational{Int}((-1)^(i)*sgn)
            end
        end
        
    
        # positive face
        key = fc[i][2]
        if key in keys(cdict)
            indSgn=cdict[key]
            ind = indSgn[1]
            sgn = indSgn[2]

            pos = findfirst(x->x==ind,c_rows)
            if pos == nothing 
                push!(c_rows,ind)
                push!(c_vals, Rational{Int}((-1)^(i+1)*sgn))
            else
                c_vals[pos]+=Rational{Int}((-1)^(i+1)*sgn)
            end
        end
    end

    return [c_rows,c_vals]
end

#= Function to calculate image of boundary map del_n: L_n -> L_n-1
    Inputs: 
        1) the array of n dimensional maps (EQCs of maps)
        2) the array of n-1 dimensional maps (EQCs of maps)
        3) the n-dimensional face list

    Outputs: 
        1) an array of sparse coordinate vectors 
    
    Notes: 
        
=#
function calculateImage(nonDegens,lowerNonDegen,faceList)
    B = Channel{}(length(nonDegens))
    cdict = coordDict(lowerNonDegen)

    @threads for map in nonDegens
        im = boundarySum(map[1],cdict,faceList) # take first rep of classes
        put!(B,im)        
    end 
    close(B)
    image = collect(B)
    return image
end

#= function to generate the matrix representing the n-th boundary map, stored sparsly
    Inputs: 
        1) the array of n dimensional maps (EQCs of maps)
        2) the array of n-1 dimensional maps (EQCs of maps)
        3) the n-dimensional face list

    Outputs: 
        1) the matrix of the n-th boundary map stored as a sparse matrix
    
    Notes: 
        
=#
function boundaryMapMatrixSparse(nonDegens,lowerNonDegens, faceList)  # for now cant do sparse, fix
    img = calculateImage(nonDegens, lowerNonDegens,faceList)

    if length(img)==0
        return []
    end

    row_inds = []
    col_inds = []
    vals = Rational{Int}[]

    for (i,im) in enumerate(img)
        for j = 1:length(im[1])
            push!(row_inds, im[1][j])
            push!(vals, Rational{Int}(im[2][j]))
            push!(col_inds,i)
        end
    end

    lenVec = length(lowerNonDegens)

    m=sparse(row_inds, col_inds, vals, lenVec, length(img))
    return m
end

#= Function to generate the matrix of the n+1 boundary map directly without storing any intermediate variables
    Inputs: 
        1) an array of n dimensional EQCs
        2) an array of n dimensional non degenerate EQCs
        3) n+1 face list
        4) the nhood dictionary of the graph
        5) the n+1 dimensional EQClist

    Outputs: 
        1) an array of sparse coordinate vectors 
    
    Notes: 
        - right now using the conjecture, if false remove the break line 
        
=#
function graphMapsMatrix(A,lowerNonDegen,facesList,nhoodDict,EQClist) # A = array of the n-1 cubes, G = graph
    
    B = Channel{}(length(A)*length(A)*length(A[1])*length(A[1])) # channel of coordinate vectors

    cdict = coordDict(lowerNonDegen)

    lenVec = length(lowerNonDegen)

    @threads for x = 1:length(A)-1  
 
        for y = x+1:length(A)
            
            h=A[x]
            k=A[y]
            
            f=h[1]
            
                for g in k
                    if isPairNcube(f,g,nhoodDict) # check if forms an n+1 cube
                        
                        # check first to see if degenerate
                        degen=checkSemiDegen(vcat(f[1],g[1]),EQClist)

                        c_rows = []
                        c_vals = Rational{Int}[]
                        
                        if !(degen) 
                            
                            fc = faces(vcat(f[1],g[1]),facesList) # generate the faces
                            
                            for i = 1:length(fc) # loop through the faces
                                
                                # negative face
                                key = fc[i][1]
                                if key in keys(cdict)
                                    indSgn=cdict[key]
                                    ind = indSgn[1]
                                    sgn = indSgn[2]

                                    pos = findfirst(x->x==ind, c_rows)

                                    if pos == nothing
                                        push!(c_rows,ind)
                                        push!(c_vals, Rational{Int}((-1)^(i)*sgn))
                                    else
                                        c_vals[pos] += Rational{Int}((-1)^(i)*sgn)
                                    end
                                end
                                
                            
                                # positive face
                                key = fc[i][2]
                                if key in keys(cdict)
                                    indSgn=cdict[key]
                                    ind = indSgn[1]
                                    sgn = indSgn[2]

                                    pos = findfirst(x->x==ind, c_rows)

                                    if pos == nothing
                                        push!(c_rows,ind)
                                        push!(c_vals, Rational{Int}((-1)^(i+1)*sgn))
                                    else
                                        c_vals[pos] += Rational{Int}((-1)^(i+1)*sgn)
                                    end
                                end
                                
                            end
                             
                            put!(B,[c_rows,c_vals])
                             
                        end
                        
                    end
                end
        end 
        
    end
    
    close(B)
    img = collect(B)
  
    if length(img)==0
        return []
    end

    row_inds = []
    col_inds = []
    vals = Rational{Int}[]

    for (i,im) in enumerate(img)
        for j = 1:length(im[1])
            push!(row_inds, im[1][j])
            push!(vals, Rational{Int}(im[2][j]))
            push!(col_inds,i)
        end
    end
    
    m=sparse(row_inds, col_inds, vals, lenVec, length(img))

    return m
end

#= Function to generate the n-th homology of a graph
    Inputs: 
        1) the graph G
        2) the dimension to compute
    Outputs: 
        1) the n-th homology group
    
    Notes: 
        
=#
function cubicalHomology(G::graph,n::Int)
    
    #Error handling
    #verify graph has number entries
    if !all(x -> x isa Number, G.vertices)
        error("Error: all vertices of the graph must be numbers. Try running the relabel_vertices function on the graph before computing homology") 
    end

    #verify positive dimension
    if n < 0
        error("Error: dimension must be at least 0") 
    end


    #computing the homology
    #Initialize dicts of EQC and faces
    EQCdict = Dict{Int, Array{Any}}()

    faceDict = Dict{Int, Array{Any}}()

    #fill the dictionaries
    for i = 1:n+1
        EQCdict[i] = generateEQClist(i)
    end

    for i = 0:n+1
        faceDict[i] = faceList(i)
    end
    
    # initialize nhood dict for checking maps
    nhoodDict = get_nhood_dict(G)

    # Initialize map the dictionary
    maps = Dict{Int, Array{Any}}()  

    # Populate the dictionary
    maps[-1]= [[[nothing,1]]] # for 0-th homology groups
    maps[0] = graphMaps([],G,[],[])
    for i=1:n

        mi=graphMaps(maps[i-1],G,EQCdict[i],nhoodDict)
        
        maps[i]=mi
    
        maps[i]=remove_duplicates(maps[i])
        
    end
    
    # remove degeneracies and semi degeneracies
    
    nonDegen = semiNonDegen(maps[n])

    lowerNonDegen = semiNonDegen(maps[n-1])

    
    # n-th boundary map 
    delN = boundaryMapMatrixSparse(nonDegen, lowerNonDegen, faceDict[n])
    
    # determine the matrix for the n+1-th boundary map
    delN1 = graphMapsMatrix(maps[n],nonDegen,faceDict[n+1],nhoodDict,EQCdict[n+1])
   
    if delN1 == [] || iszero(delN1) # zero matrix becomes wrong type and can cause issues
        dimIM = 0
    else
        dimIM=rank(delN1) # faster to do transpose
    end
    
    
    # handling for different cases
    if delN==[]
        dimKer = 0
    elseif iszero(delN) # zero matrix becomes wrong type and can cause issues
        dimKer = size(delN,2)
    else
        dimKer=size(delN,2)-rank(delN) # faster to do transpose
    end


    # homology dimension
    return dimKer-dimIM
    
end

#=================================================================================
# functions for fast H0 and H1 computations over Z2 
=================================================================================#

# gets the corresponding vector of the boundary of map
function get_vec(map, cdict, len)
    v = zeros(ZZ2,len)
    faces = [[map[1],map[2]],[map[3],map[4]],[map[2],map[4]],[map[1],map[3]]]

    for i = 1:4
        try
            c = cdict[faces[i]]
            v[c] += 1
        catch e 
        end
    end

    return v
end

# returns H0, H1
function computeH1(G) 

    G = relabel_vertices(G)
    # initialize nhood dict for checking maps
    nhoodDict = get_nhood_dict(G)

    C0 = G.vertices
    C1 = [[e[1], e[2]] for e in G.edges if e[1] != e[2]]
    C1 = unique!(C1)

    # n-th boundary map 
    B = Channel{}(length(C1))
    len = length(G.vertices)

    @threads for c in C1 
        v = zeros(ZZ2,len)
        v[[c[1]]] .+= 1
        v[[c[2]]] .+= 1
        put!(B,v)
    end

    close(B)
    img = collect(B)
  
    if length(img)==0
        delN = []
    else
        delN=hcat(img ...)
    end
    
    
    # determine the 2 boundary matrix
    B = Channel{}(length(G.edges)*length(G.edges)) # channel of coordinate vectors

    cdict = Dict{}()
    
    for i = 1:length(C1)
        cdict[C1[i]] = i
    end

    len = length(C1)

    @threads for i = 1:length(G.vertices)
        v = G.vertices[i]
        nv = nhoodDict[v]

        geq_v = G.vertices[i+1:end]
        
        #Handling for all distinct vertices
        for w in intersect(nv, geq_v)
            nw = nhoodDict[w]
            j = findfirst(x -> x == w, G.vertices)
            geq_w = G.vertices[j+1:end]
            for v2 in intersect(nv, geq_w)
                nv2 = nhoodDict[v2] 
                for w2 in intersect(nw, nv2, geq_v)
                    
                    vec = get_vec([v,w,v2,w2],cdict,len)
                    put!(B,vec)
                    vec = get_vec([v,v2,w,w2],cdict,len)
                    put!(B,vec)
    
                    
                end
                
                # if w2 = v 
                vec = get_vec([v,w,v2,v],cdict,len)
                put!(B,vec)

                vec = get_vec([v,v2,w,v],cdict,len)
                put!(B,vec)

            end

            # if v2 = v 
            for w2 in intersect(nw, nv, geq_w)

                vec = get_vec([v,w,v,w2],cdict,len)
                put!(B,vec)

                vec = get_vec([v,v,w,w2],cdict,len)
                put!(B,vec)
    

            end
            
            for w2 in filter!(x -> x != w, intersect(nw, geq_v)) 
                vec = get_vec([w,w2,v,w],cdict,len)
                put!(B,vec)
                vec = get_vec([w,v,w2,w],cdict,len)
                put!(B,vec)
            end

            vec = get_vec([v,v,w,v],cdict,len)
            put!(B,vec)

            vec = get_vec([w,w,v,w],cdict,len)
            put!(B,vec)

        end


    end
    
    close(B)
    img = collect(B)
  
    delN1=hcat(img ...)


    if delN1 == []
        dimIM = 0
    else
        dimIM = rank!(delN1)
    end

    if delN == []
        dimIM1 = 0
    else
        dimIM1 = rank!(delN)
    end


    # handling for different cases
    if delN==[]
        dimKer = 0
    else
        dimKer=size(delN,2)-dimIM1
    end


    # homology dimension H0, H1
    return length(G.vertices) - dimIM1, dimKer-dimIM
    
end


#===================================================================
# Functions to compute the homology of a Graph mod p where p is the smallest prime larger than n+1
===================================================================#

#= Function computes the boundary of a map and stores it as a sparse vector
    Inputs: 
        1) the map as an array
        2) the lower non degenerate cubes
        3) the n-dimensional face list

    Outputs: 
        1) the boundary as a coordinate vector, stored sparsly
    
    Notes: 
        
=#
function boundarySumModP(map,cdict,faceList,p)
    
    c_rows = []
    c_vals = Mod{p}[]
    
    fc = faces(map[1],faceList) # generate faces

    for i = 1:length(fc) # loop through faces
        # negative face
        key = fc[i][1]
        if key in keys(cdict)
            indSgn=cdict[key]
            ind = indSgn[1]
            sgn = indSgn[2]

            pos = findfirst(x->x==ind,c_rows)
            if pos == nothing 
                push!(c_rows,ind)
                push!(c_vals,Mod{p}((-1)^(i)*sgn))
            else
                c_vals[pos] += Mod{p}((-1)^(i)*sgn)
            end
        end
        
    
        # positive face
        key = fc[i][2]
        if key in keys(cdict)
            indSgn=cdict[key]
            ind = indSgn[1]
            sgn = indSgn[2]

            pos = findfirst(x->x==ind,c_rows)
            if pos == nothing 
                push!(c_rows,ind)
                push!(c_vals,Mod{p}((-1)^(i+1)*sgn))
            else
                c_vals[pos] += Mod{p}((-1)^(i+1)*sgn)
            end
        end
    end

    return [c_rows,c_vals]
end

#= Function to calculate image of boundary map del_n: L_n -> L_n-1
    Inputs: 
        1) the array of n dimensional maps (EQCs of maps)
        2) the array of n-1 dimensional maps (EQCs of maps)
        3) the n-dimensional face list

    Outputs: 
        1) an array of sparse coordinate vectors 
    
    Notes: 
        
=#
function calculateImageModP(nonDegens,lowerNonDegen,faceList,p)
    B = Channel{}(length(nonDegens))
    cdict = coordDict(lowerNonDegen)

    @threads for map in nonDegens
        im = boundarySumModP(map[1],cdict,faceList,p) # take first rep of classes
        put!(B,im)        
    end 
    close(B)
    image = collect(B)
    return image
end

#= function to generate the matrix representing the n-th boundary map, stored sparsly
    Inputs: 
        1) the array of n dimensional maps (EQCs of maps)
        2) the array of n-1 dimensional maps (EQCs of maps)
        3) the n-dimensional face list

    Outputs: 
        1) the matrix of the n-th boundary map stored as a sparse matrix
    
    Notes: 
        
=#
function boundaryMapMatrixSparseModP(nonDegens,lowerNonDegens, faceList,p)  # for now cant do sparse, fix
    img = calculateImageModP(nonDegens, lowerNonDegens,faceList,p)

    if length(img)==0
        return []
    end

    row_inds = []
    col_inds = []
    vals = Mod{p}[]

    for (i,im) in enumerate(img)
        for j = 1:length(im[1])
            push!(row_inds, im[1][j])
            push!(vals, Mod{p}(im[2][j]))
            push!(col_inds,i)
        end
    end
    
    lenVec = length(lowerNonDegens)

    m=sparse(row_inds, col_inds, vals, lenVec, length(img))
    
    return m
end

#= Function to generate the matrix of the n+1 boundary map directly without storing any intermediate variables
    Inputs: 
        1) an array of n dimensional EQCs
        2) an array of n dimensional non degenerate EQCs
        3) n+1 face list
        4) the nhood dictionary of the graph
        5) the n+1 dimensional EQClist

    Outputs: 
        1) an array of sparse coordinate vectors 
    
    Notes: 
        - right now using the conjecture, if false remove the break line 
        
=#
function graphMapsMatrixModP(A,lowerNonDegen,facesList,nhoodDict,EQClist,p) # A = array of the n-1 cubes, G = graph
    
    B = Channel{}(length(A)*length(A)*length(A[1])*length(A[1])) # channel of coordinate vectors

    cdict = coordDict(lowerNonDegen)

    lenVec = length(lowerNonDegen)

    @threads for x = 1:length(A)-1  
 
        for y = x+1:length(A)
            
            h=A[x]
            k=A[y]
            
            f=h[1]
            
                for g in k
                    if isPairNcube(f,g,nhoodDict) # check if forms an n+1 cube
                        
                        # check first to see if degenerate
                        degen=checkSemiDegen(vcat(f[1],g[1]),EQClist)
                        
                        if !(degen) 
                            c_rows = []
                            c_vals = Mod{p}[]
                            
                            fc = faces(vcat(f[1],g[1]),facesList) # generate the faces
                            
                            for i = 1:length(fc) # loop through the faces
                                
                                # negative face
                                key = fc[i][1]
                                if key in keys(cdict)
                                    indSgn=cdict[key]
                                    ind = indSgn[1]
                                    sgn = indSgn[2]

                                    pos = findfirst(x->x==ind, c_rows)
                                    if pos == nothing
                                        push!(c_rows,ind)
                                        push!(c_vals, Mod{p}((-1)^(i)*sgn))
                                    else
                                        c_vals[pos] +=  Mod{p}((-1)^(i)*sgn)
                                    end
                                end
                                
                            
                                # positive face
                                key = fc[i][2]
                                if key in keys(cdict)
                                    indSgn=cdict[key]
                                    ind = indSgn[1]
                                    sgn = indSgn[2]

                                    pos = findfirst(x->x==ind, c_rows)
                                    if pos == nothing
                                        push!(c_rows,ind)
                                        push!(c_vals, Mod{p}((-1)^(i+1)*sgn))
                                    else
                                        c_vals[pos] +=  Mod{p}((-1)^(i+1)*sgn)
                                    end
                                end
                                
                            end
                             
                            put!(B,[c_rows,c_vals])
                             
                        end
                        
                        
                    end
                end
            
        end 
        
    end
    
    close(B)
    img = collect(B)
  
    if length(img)==0
        return []
    end
    
    row_inds = []
    col_inds = []
    vals = Mod{p}[]

    for (i,im) in enumerate(img)
        for j = 1:length(im[1])
            push!(row_inds, im[1][j])
            push!(vals, Mod{p}(im[2][j]))
            push!(col_inds,i)
        end
    end
    

    m=sparse(row_inds, col_inds, vals, lenVec, length(img))

    return m
end

#= Function to generate the n-th homology of a graph
    Inputs: 
        1) the graph G
        2) the dimension to compute
    Outputs: 
        1) the n-th homology group
    
    Notes: 
        
=#
function cubicalHomologyModP(G::graph,n::Int)
    
    #Error handling
    #verify graph has number entries
    if !all(x -> x isa Number, G.vertices)
        error("Error: all vertices of the graph must be numbers. Try running the relabel_vertices function on the graph before computing homology") 
    end

    #verify positive dimension
    if n < 0
        error("Error: dimension must be at least 0") 
    end

    # determine which field to compute over
    p = nextprime(n+2)


    #computing the homology
    #Initialize dicts of EQC and faces
    EQCdict = Dict{Int, Array{Any}}()

    faceDict = Dict{Int, Array{Any}}()

    #fill the dictionaries
    for i = 1:n+1
        EQCdict[i] = generateEQClist(i)
    end

    for i = 0:n+1
        faceDict[i] = faceList(i)
    end
    
    # initialize nhood dict for checking maps
    nhoodDict = get_nhood_dict(G)

    # Initialize map the dictionary
    maps = Dict{Int, Array{Any}}()  

    # Populate the dictionary
    maps[-1]= [[[nothing,1]]] # for 0-th homology groups
    maps[0] = graphMaps([],G,[],[])
    for i=1:n

        mi=graphMaps(maps[i-1],G,EQCdict[i],nhoodDict)
        
        maps[i]=mi
    
        maps[i]=remove_duplicates(maps[i])
        
    end
    
    # remove degeneracies and semi degeneracies
    
    nonDegen = semiNonDegen(maps[n])

    lowerNonDegen = semiNonDegen(maps[n-1])

    
    # n-th boundary map 
    delN = boundaryMapMatrixSparseModP(nonDegen, lowerNonDegen, faceDict[n],p)
    
    # determine the matrix for the n+1-th boundary map
    delN1 = graphMapsMatrixModP(maps[n],nonDegen,faceDict[n+1],nhoodDict,EQCdict[n+1],p)


    if delN1 == [] || iszero(delN1) # zero matrix becomes wrong type and can cause issues
        dimIM = 0
    else
        dimIM=ref_rank!(delN1,true) # faster to do transpose
    end
    
    
    # handling for different cases
    if delN==[]
        dimKer = 0
    elseif iszero(delN) # zero matrix becomes wrong type and can cause issues
        dimKer = size(delN,2)
    else
        dimKer=size(delN,2)-ref_rank!(delN,true) # faster to do transpose
    end


    # homology dimension
    return dimKer-dimIM
    
end

#===================================================================#
# JSON File Reading and Benchmarking Functions
#===================================================================#

# Convert adjacency list to graph structure
function adjacency_to_graph(vertices, adj)
    edges = []
    
    # Convert adjacency list to edge list
    for (v, neighbors) in adj
        for neighbor in neighbors
            push!(edges, (v, neighbor))
        end
    end
    
    # Make symmetric and reflexive
    edges = makeSymmetricReflexive(edges)
    
    return graph(vertices, edges)
end

# Calculate H0 with memory management for cubical homology
function calculate_cubical_h0(vertices, adj; mod_p=false, dimension=0)
    try
        # Convert to graph structure
        G = adjacency_to_graph(vertices, adj)
        G = relabel_vertices(G)
        
        result = if !mod_p
            cubicalHomology(G, dimension)
        else
            computeH1(G)[1]
        end
        
        # Force garbage collection to free memory
        G = nothing
        GC.gc()
        
        return result
    catch e
        if isa(e, OutOfMemoryError)
            println("Out of memory error in cubical H$dimension calculation")
            return "OOM"
        else
            rethrow(e)
        end
    end
end

# Calculate H1 via edge graph
function calculate_cubical_h1(vertices, adj; mod_p=false)
    try
        # Convert to graph structure
        G = adjacency_to_graph(vertices, adj)
        G = relabel_vertices(G)
        
        result = if !mod_p
            # Computes rational (cubical) 1st homology
            cubicalHomology(G, 1)
        else
            # Computes 1st homology via the edge graph method.
            computeH1(G)[2]
        end
        
        # Force garbage collection to free memory
        G = nothing
        GC.gc()
        
        return result
    catch e
        if isa(e, OutOfMemoryError)
            println("Out of memory error in cubical H1 calculation")
            return "OOM"
        else
            rethrow(e)
        end
    end
end

# Get the number of graphs without loading all data
function get_graph_count(filename="graphs.json")
    base_dir = dirname(@__FILE__)
    path = joinpath(base_dir, "..", "data", filename)
    
    # Read and parse just to count graphs
    data = JSON3.read(read(path, String))
    count = length(data)
    
    # Clear the data immediately
    data = nothing
    GC.gc()
    
    return count
end

# Load a single graph by index
function load_single_graph(filename, index)
    base_dir = dirname(@__FILE__)
    path = joinpath(base_dir, "..", "data", filename)
    
    # Read the file content
    content = read(path, String)
    data = JSON3.read(content)
    
    # Extract the specific graph
    if index <= length(data)
        graph = data[index]
        
        # Clear everything else from memory immediately
        data = nothing
        content = nothing
        GC.gc()
        
        return graph
    else
        data = nothing
        content = nothing
        GC.gc()
        return nothing
    end
end

# Check if a graph should be processed (size filtering)
function should_process_graph(graph, max_vertices=1000)
    vertices = graph["vertices"]
    return length(vertices) <= max_vertices
end

function run_cubical_homology()
    # Field selection
    println("Compute cubical homology over Q or Z/p? (q/p): ")
    input_str = "p"
    mod_p = input_str == "p"
    
    println("Number of threads available: ", nthreads())
    
    filename = "std_var_edge_high_verts_4.json"
    total_graphs = get_graph_count(filename)
    println("Total graphs in file: $total_graphs")
    
    processed_count = 0
    skipped_count = 0
    results = []
    
    # Process graphs one at a time
    for i in 1:total_graphs
        # Load only this graph
        graph = load_single_graph(filename, i)
        
        if graph === nothing
            println("Error loading graph $i")
            continue
        end
        
        # Check if we should process this graph
        if !should_process_graph(graph)
            skipped_count += 1
            graph = nothing
            GC.gc()
            continue
        end
        
        processed_count += 1
        
        name = get(graph, "name", "Graph_$i")
        vertices = graph["vertices"]
        
        # Convert everything to strings
        vertices_str = String.(vertices)
        adj_raw = Dict(k => v for (k, v) in graph["adjacency_list"])
        adj = Dict(String(k) => String.(v) for (k, v) in adj_raw)
        
        println("Processing graph $processed_count (index $i): $name ($(length(vertices_str)) vertices)")
        
        
        # Clear intermediate data and force garbage collection
        vertices_str = nothing
        adj = nothing
        adj_raw = nothing
        vertices = nothing
        GC.gc()
        
        # Reload graph for H1 calculation (since we cleared the data)
        graph = load_single_graph(filename, i)
        vertices = graph["vertices"]
        vertices_str = String.(vertices)
        adj_raw = Dict(k => v for (k, v) in graph["adjacency_list"])
        adj = Dict(String(k) => String.(v) for (k, v) in adj_raw)
        
        # Calculate H1 with timing
        start_time_h1 = time()
        H1 = calculate_cubical_h1(vertices_str, adj, mod_p=mod_p)
        elapsed_h1 = round(time() - start_time_h1, digits=15)
        
        total_elapsed = round(elapsed_h1, digits=15)
        
        # Store results
        push!(results, (name, H1, elapsed_h1, total_elapsed))
        
        # Clear all graph data from memory
        graph = nothing
        vertices = nothing
        vertices_str = nothing
        adj = nothing
        adj_raw = nothing
        GC.gc()
        
        # Print progress
        println("  H1: $H1 ($(elapsed_h1)s)")
        println("  Total: $(total_elapsed)s")
        println("  Memory cleaned up\n")
        
        # Small delay to help with memory management
        sleep(0.1)
    end
    
    println("Processed: $processed_count graphs")
    println("Skipped (too large): $skipped_count graphs")
    
    # Print final summary
    println("="^50)
    println("FINAL RESULTS:")
    println("="^50)
    for (name, H1, elapsed_h1, total_elapsed) in results
        println("Graph: $name")
        println("H1: $H1")
        println("H1 Time: $elapsed_h1 sec")
        println("Total Time: $total_elapsed sec\n")
    end
end

function create_T3()
    vert=[0,1,2,3,4,5]
    edge=[(0,1),(1,2),(2,3),(3,4),(4,5)]
    edge=makeSymmetricReflexive(edge)
    I5=graph(vert,edge)

    Graph=multBoxProd([I5,I5,I5])

    for i=0:5
        for j=0:5
            Graph=quotient_vertices(Graph,(0,i,j),(5,i,j))
        end 
    end

    for i=0:4
        for j=0:5
            Graph=quotient_vertices(Graph,(i,0,j),(i,5,j))
        end 
    end

    for i=0:4
        for j=0:4
            Graph=quotient_vertices(Graph,(i,j,0),(i,j,5))
        end
    end
    
    return relabel_vertices(Graph)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_cubical_homology()
end