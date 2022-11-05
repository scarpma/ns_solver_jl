module IO
using Printf
using WriteVTK

using WriteVTK
using Printf

function writeSurface(
        filename::String,
        cooOfVert::Union{Matrix{Float64}, Matrix{Float32}},
        vertOfTria::Matrix{Int64},
        fieldNames::Vector{String}=Vector{String}(undef,0),
        fields...)
    if vertOfTria == nothing
        cells = [MeshCell(VTKCellTypes.VTK_VERTEX, (i, )) for i = 1:size(cooOfVert,2)]
    else
        cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE,
                    (vertOfTria[1, i], vertOfTria[2, i], vertOfTria[3, i]))
                    for i = 1:size(vertOfTria,2)
                ]
    end
        mesh = vtk_grid(filename, cooOfVert[1,:], cooOfVert[2,:], cooOfVert[3,:], cells)
    if size(fields,1) > 0
        @assert size(fieldNames,1) == size(fields,1)
        for i in eachindex(fields)
            mesh[fieldNames[i]] = fields[i]
        end
    end
    vtk_save(mesh)
    return
end

function read_gts(filename::String, T::DataType)
    lines = readlines(filename)
    nNodes, nEdges, nTria = parse.(Int64, split(lines[1], " "))
    coords = reshape(parse.(T, split(join(lines[2:nNodes+1], " "), " ")), (3,nNodes))
    vertOfEdge = reshape(parse.(Int64, split(join(lines[nNodes+2:nNodes+2+nEdges-1], " "), " ")), (2,nEdges))
    edgeOfTria = reshape(parse.(Int64, split(join(lines[nNodes+2+nEdges:end], " "), " ")), (3,nTria))

    vertOfTria = zeros(Int64, (3,nTria))
    for i=1:nTria
        e1, e2, e3 = edgeOfTria[:,i] # e3 is not used
        if vertOfEdge[2,e1] == vertOfEdge[1,e2]      #-----#-----#-----#
           v1 = vertOfEdge[1,e1] #                   1 e1  21 e2 2? e3 ?
           v2 = vertOfEdge[2,e1] 
           v3 = vertOfEdge[2,e2]
        elseif vertOfEdge[2,e1] == vertOfEdge[2,e2]  #-----#-----#-----#
           v1 = vertOfEdge[1,e1] #                   1 e1  22 e2 1? e3 ? 
           v2 = vertOfEdge[2,e1] 
           v3 = vertOfEdge[1,e2]
        elseif vertOfEdge[1,e1] == vertOfEdge[1,e2]  #-----#-----#-----#
           v1 = vertOfEdge[2,e1] #                   2 e1  11 e2 2? e3 ?
           v2 = vertOfEdge[1,e1] 
           v3 = vertOfEdge[2,e2]
        else                                         #-----#-----#-----#
           v1 = vertOfEdge[2,e1] #                   2 e1  12 e2 1? e3 ?
           v2 = vertOfEdge[1,e1] 
           v3 = vertOfEdge[1,e2]
        end
        vertOfTria[:,i] .= v1, v2, v3
    end
    return coords, vertOfTria, vertOfEdge, edgeOfTria
end

function read_gts_plus(filename::String, T::DataType)
    coords, vertOfTria, vertOfEdge, edgeOfTria = read_gts(filename, T)
    nNodes = size(coords, 2)
    nEdges = size(vertOfEdge, 2)
    nTria = size(vertOfTria, 2)
    
    nEdgeVert = zeros(Int64, (nNodes))
    midPointEdge = zeros(T, (3, nEdges))
    midPointEdge = zeros(T, (3, nEdges))
    dirEdge = zeros(T, (3, nEdges))
    lenEdge = zeros(T, nEdges)
    for i=1:nEdges
        v1, v2 = vertOfEdge[:,i]
        nEdgeVert[v1] += 1
        nEdgeVert[v2] += 1
        midPointEdge[:,i] .= 0.5.*(coords[:,v1] .+ coords[:,v2])
        dirEdge[:,i] .= (coords[:,v1] .- coords[:,v2])
        lenEdge[i] = sqrt(sum(dirEdge[:,i].^2))
        dirEdge[:,i] .= dirEdge[:,i] ./ lenEdge[i]
    end
    
    maxNEdgeVert = maximum(nEdgeVert)
    vertOfVert = zeros(Int64, (maxNEdgeVert, nNodes))
    edgeOfVert = zeros(Int64, (maxNEdgeVert, nNodes))
    nEdgeVert .= 0
    for i=1:nEdges
        v1, v2 = vertOfEdge[:,i]
        nEdgeVert[v1] += 1
        nEdgeVert[v2] += 1
        vertOfVert[nEdgeVert[v1], v1] = v2
        vertOfVert[nEdgeVert[v2], v2] = v1
        edgeOfVert[nEdgeVert[v1], v1] = i 
        edgeOfVert[nEdgeVert[v2], v2] = i
    end
    return coords, vertOfTria, vertOfEdge, vertOfVert,
    edgeOfTria, nEdgeVert, midPointEdge, lenEdge, dirEdge
end

function test_read_gts()
    coords, vertOfTria, vertOfEdge,
    vertOfVert, edgeOfTria, nEdgeVert,
    midPointEdge, lenEdge, dirEdge = read_gts_plus("./aorta1.gts", Float32);
    coords, vertOfTria, = read_gts("./aorta1.gts", Float32);
    writeSurface("prova", coords, vertOfTria, nEdgeVert)
    writeSurface("midPoints", midPointEdge, nothing, dirEdge, lenEdge)
end

function staggeredToNormal(staggered)
    n1, n2 = size(staggered)
    n1 = n1 - 1
    n2 = n2 - 1
    normal = Array{eltype(staggered)}(undef, (n1,n2))
    for i=1:n1
        for j=1:n2
            normal[i,j] = 1/4*(staggered[i,j]+staggered[i,j+1]+staggered[i+1,j]+staggered[i+1,j+1])
        end
    end
    return normal
end

function writeUVPFiedsToFile(filename::String,x1,x2,u,v,p)
    h1 = x1[2] - x1[1]
    n1 = length(x1)
    x1n = collect(range(-h1/2, (n1+0.5)*h1, n1+1))
    println(x1n)
    h2 = x2[2] - x2[1]
    n2 = length(x2)
    x2n = collect(range(-h2/2, (n2+0.5)*h2, n2+1))
    vtkfile = vtk_grid(filename, x1n, x2n)
    vtkfile["U"] = vcat(u,v)
    vtkfile["p"] = staggeredToNormal(p)
    outfiles = vtk_save(vtkfile)
end

function writeSingleFieldToFile(filename::String,x1,x2,p,fieldName::String)
    h1 = x1[2] - x1[1]
    n1 = length(x1)
    x1n = collect(range(-h1/2, (n1+0.5)*h1, n1+1))
    h2 = x2[2] - x2[1]
    n2 = length(x2)
    x2n = collect(range(-h2/2, (n2+0.5)*h2, n2+1))
    vtkfile = vtk_grid(filename, x1n, x2n)
    vtkfile[fieldName] = staggeredToNormal(p)
    outfiles = vtk_save(vtkfile)
end

function writeFieldsToFile(filename::String,x1,x2,fields...)
    # fields values are assumed to be on cells
    # => create grid points x1n and x2n among grid cells
    # x1 and x2 are on grid cells (centroids)
    h1 = x1[2] - x1[1]
    n1 = length(x1)
    x1n = collect(range(-h1/2, (n1+0.5)*h1, n1+1))
    h2 = x2[2] - x2[1]
    n2 = length(x2)
    x2n = collect(range(-h2/2, (n2+0.5)*h2, n2+1))
    vtkfile = vtk_grid(filename, x1n, x2n)
    for i in eachindex(fields)
        vtkfile[@sprintf("f%2d",i)] = fields[i]
    end
    outfiles = vtk_save(vtkfile)
end

function readSTLdim(filename::String)::Int64
    lines = readlines(filename)
    ntria::Int64 = 0
    i::Int64 = 1
    while true
        if occursin("endsolid", lines[i])
            break
        elseif occursin("outer ", lines[i])
            ntria = ntria + 1
            i += 4
        else
            i += 1
        end
    end
    return ntria
end

function readSTL(filename::String)::Array{Float64, 2}
    ntria = readSTLdim(filename)
    vertOfTria = Array{Float64, 2}(undef, (3,3*ntria))
    lines = readlines(filename)
    i::Int64 = 1
    j::Int64 = 1
    while true
        if occursin("endsolid", lines[i])
            break
        end
        if occursin("outer ", lines[i])
            line = strip(lines[i])
            vertOfTria[:,j  ] = parse.(Float64, split(strip(lines[i+1]), " ")[2:end])
            vertOfTria[:,j+1] = parse.(Float64, split(strip(lines[i+2]), " ")[2:end])
            vertOfTria[:,j+2] = parse.(Float64, split(strip(lines[i+3]), " ")[2:end])
            j = j + 3
            i += 4
            else
                i += 1
            end
        end
    return vertOfTria
end

function stlConnectivity(vertOfTria::Array{Float64, 2})
    """
    this function iterates over vertices found inside an STL
    file and keeps only new vertices inside the array <nodes>.
    Each new vertex not closer than <tol> to any other vertices
    in <nodes> is added to <nodes>.
    
    <whois> maps from unordered vertices to unordered vertices
        all_vertices_index --> corresponding unique vertex inside all_vertices_index
    <nodes> maps from ordered vertices to unordered vertices
        unique vertex index from ordered vertex array --> all_vertices_index
    <nodes_inv> maps from unordered vertices to ordered vertices
        all_vertices_index --> corresponding unique vertex inside ordered vertex array
    
    """    
    tol = 1.0e-12 # m^2
    ntria = div(size(vertOfTria, 2), 3)
    maxNOfVert = 3*ntria # temp dim for new ordered vertices array
    alreadyAssigned = zeros(Int64, (ntria*3))
    alreadyAssigned[1] = 1
    whois = Array{Int64, 1}(undef, (3*ntria))
    nodes = Array{Int64, 1}(undef, (maxNOfVert))
    nodes_inv = Array{Int64, 1}(undef, (3*ntria))
    nodes[1] = 1
    nodes_inv[1] = 1
    whois[1] = 1
    whois[1] = 1
    nVertOrdered = 1 # current len of ordered vertices <nodes>
    last_j::Int64 = 0
    for i=1:3*ntria
        if alreadyAssigned[i] == 0
            dist_sum = 0.::Float64
            isDouble = false
            for j=1:nVertOrdered
                dist2 = sum((vertOfTria[:,i] - vertOfTria[:,nodes[j]]).^2)
                dist_sum += dist2
                if dist2 < tol
                    # vertex i is closer than tol to vertex nodes[j]
                    # this is not a new vertex !
                    whois[i] = nodes[j]
                    alreadyAssigned[i] = 1
                    nodes_inv[i] = j
                    isDouble = true
                    break
                end
            end
            # vertex <i> is not closer than <tol> to any vertex inside <nodes>
            # let's add vertex <i> to <nodes>
            if isDouble == false
                nVertOrdered += 1
                nodes[nVertOrdered] = i
                nodes_inv[i] = nVertOrdered
                whois[i] = i
                alreadyAssigned[i] = 1    
            end
        end
    end
    resize!(nodes, nVertOrdered)
    return nodes, nodes_inv, whois
end

end
