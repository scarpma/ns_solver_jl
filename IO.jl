module IO
using Printf
using WriteVTK

using WriteVTK
using Printf

function writeSurface(filename, cooOfVert, vertOfTria, fields...)
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
        for i in eachindex(fields)
            mesh[@sprintf("f%2d",i)] = fields[i]
        end
    end
    vtk_save(mesh)
    return
end

function read_gts(filename, T)
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

function read_gts_plus(filename, T)
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

end
