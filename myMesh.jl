module myMesh

include("./IO.jl")
using .IO
using LinearAlgebra

export SurfaceMesh, compute_distances!, compute_angle!,
       compute_triaNormalAndArea!

mutable struct SurfaceMesh{T}
    nv::Int64 # number of vertices
    ne::Int64 # number of edges
    nt::Int64 # number of triangles
    xyzOfVert::Array{T, 2}
    vertOfTria::Array{Int64, 2}
    vertOfEdge::Array{Int64, 2}
    edgeOfTria::Array{Int64, 2}
    triaOfEdge::Array{Int64, 2}
    nTriaOfEdge::Array{Int64, 1}
    vert1234OfEdge::Array{Int64, 2}
    lenOfEdges0::Array{T, 1}
    triaNormals0::Array{T,2}
    triaAreas0::Array{T,1}
    angleOfEdges0::Array{T,1}
    function SurfaceMesh(
        xyzOfVert::Union{Array{Float64, 2}, Array{Float32, 2}},
        vertOfTria::Array{Int64, 2},
        vertOfEdge::Array{Int64, 2},
        edgeOfTria::Array{Int64, 2},
    )
        # initialize arrays
        T = eltype(xyzOfVert)
        nv = size(xyzOfVert, 2)
        ne = size(vertOfEdge, 2)
        nt = size(vertOfTria, 2)
        triaOfEdge = zeros(Int64, (2,ne))
        nTriaOfEdge = zeros(Int64, ne)
        vert1234OfEdge = zeros(Int64, (4, ne))
        lenOfEdges0 = zeros(T, ne)
        triaNormals0 = zeros(T, (3,nt))
        triaAreas0 = zeros(T, nt)
        angleOfEdges0 = zeros(T, ne)
        
        # create instance
        m = new{eltype(xyzOfVert)}(
            nv, ne, nt, xyzOfVert, vertOfTria,
            vertOfEdge, edgeOfTria, triaOfEdge, 
            nTriaOfEdge, vert1234OfEdge, lenOfEdges0,
            triaNormals0, triaAreas0, angleOfEdges0)
        
        # populate arrays
        compute_vert1234OfEdge!(m.vert1234OfEdge, m)
        compute_triaOfEdge!(m.triaOfEdge, m.nTriaOfEdge, m)
        compute_distances!(m.lenOfEdges0, m)
        compute_triaNormalAndArea!(m.triaNormals0, m.triaAreas0, m)
        compute_angle!(m.angleOfEdges0, m, m.triaNormals0)
        
        return m
    end
end

function compute_triaOfEdge!(triaOfEdge, nTriaOfEdge, m)
    for i=1:m.nt
        e1, e2, e3 = m.edgeOfTria[:,i]
        nTriaOfEdge[e1] += 1
        nTriaOfEdge[e2] += 1
        nTriaOfEdge[e3] += 1
        triaOfEdge[nTriaOfEdge[e1], e1] = i
        triaOfEdge[nTriaOfEdge[e2], e2] = i
        triaOfEdge[nTriaOfEdge[e3], e3] = i
    end
end

function compute_vert1234OfEdge!(vert1234OfEdge, m)
    for i=1:m.ne
        t1, t2 = m.triaOfEdge[:,i]
        v2, v3 = m.vertOfEdge[:,i]
        if (t1 != 0 && t2 != 0)
            v1_, v2_, v3_ = m.vertOfTria[:,t1]
            if (v1_!=v2 && v1_!=v3)
                v1 = v1_
            elseif (v2_!=v2 && v2_!=v3)
                v1 = v2_
            else
                v1 = v3_
            end
            v1_, v2_, v3_ = m.vertOfTria[:,t2]
            if (v1_!=v2 && v1_!=v3)
                v4 = v1_
            elseif (v2_!=v2 && v2_!=v3)
                v4 = v2_
            else
                v4 = v3_
            end
            vert1234OfEdge[:,i] .= v1, v2, v3, v4
        end
    end
end

function compute_triaOfVert(m)
    nTriaOfVert = zeros(Int64, m.nv)
    for i=1:m.nt
        v1, v2, v3 = m.vertOfTria[:,i]
        nTriaOfVert[v1] += 1
        nTriaOfVert[v2] += 1
        nTriaOfVert[v3] += 1
    end
    triaOfVert = zeros(Int64, (maximum(nTriaOfVert), m.nv))
    nTriaOfVert = zeros(Int64, m.nv)
    for i=1:m.nt
        v1, v2, v3 = m.vertOfTria[:,i]
        nTriaOfVert[v1] += 1
        nTriaOfVert[v2] += 1
        nTriaOfVert[v3] += 1
        triaOfVert[nTriaOfVert[v1], v1] = i
        triaOfVert[nTriaOfVert[v2], v2] = i
        triaOfVert[nTriaOfVert[v3], v3] = i
    end
    return triaOfVert, nTriaOfVert
end


function compute_distances!(dist, m)
    for i=1:m.ne
        v1, v2 = m.vertOfEdge[:,i]
        dist[i] = sqrt(sum((m.xyzOfVert[:,v1] - m.xyzOfVert[:,v2]).^2.))
    end
    return dist
end

function compute_triaNormalAndArea!(normal, area, m)
    for i=1:m.nt
        v1 = m.xyzOfVert[:,m.vertOfTria[1,i]]
        v2 = m.xyzOfVert[:,m.vertOfTria[2,i]]
        v3 = m.xyzOfVert[:,m.vertOfTria[3,i]]
        v21 = v2 .- v1
        v31 = v3 .- v1
        normal[:,i] = LinearAlgebra.cross(v21, v31)
        area[i] = sqrt(sum(normal[:,i].^2.)) / 2.
        normal[:,i] = normal[:,i] ./ (area[i]*2)
    end
end

function cell_array_to_point_array(arr_, m)
    if ndims(arr_) == 1
        arr = reshape(arr_, (1,size(arr_,1)))
    else
        arr = arr_
    end
    triaOfVert, nTriaOfVert = compute_triaOfVert(m)
    normals = zeros((3,m.nt))
    areas = zeros(m.nt)
    areaOfVert = zeros(m.nv)
    compute_triaNormalAndArea!(normals, areas, m)
    parr = zeros(eltype(arr), (size(arr,1), m.nv))
    for i=1:m.nv
        for j=1:nTriaOfVert[i]
            t = triaOfVert[j,i]
            parr[:,i] += arr[:,t] .* areas[t]
            areaOfVert[i] += areas[t]
        end
    end
    for i=1:m.nv
        parr[:,i] = parr[:,i] ./ areaOfVert[i]
    end
    return parr
end

function compute_angle!(theta, m::SurfaceMesh, normals)
    for i=1:m.ne
        t1 = m.triaOfEdge[1,i]
        t2 = m.triaOfEdge[2,i]
        if (t1!=0 && t2!=0)
            n1 = normals[:,t1]
            n2 = normals[:,t2]
            pv = cross(n1, n2)
            tmp1 = sqrt(sum(pv.^2.))
            tmp2 = sum(n1.*n2)
            theta[i]= atan(tmp1,tmp2)
        end
    end
end

function edgeArrayToTriaArray(m::SurfaceMesh, arr)
    T = eltype(m.xyzOfVert)
    newArr = zeros(T, m.nt)
    for i=1:m.ne
        t1, t2 = m.triaOfEdge[:,i]
        newArr[t1] += arr[i]
        if m.nTriaOfEdge[i] > 1
            newArr[t2] += arr[i]
        end
    end
    newArr = newArr ./ 3.
    return newArr
end

function test_compute_angle()
    coords, vertOfTria, vertOfEdge,
    vertOfVert, edgeOfTria, nEdgeVert, 
    midPointEdge, lenEdge, dirEdge = IO.read_gts_plus("./aorta1.gts", Float64)

    m = SurfaceMesh(coords, vertOfTria, vertOfEdge, edgeOfTria)
    thetasDegTria = edgeArrayToTriaArray(m, m.angleOfEdges0) .* 360. / 2pi
    IO.writeSurface("provatheta", m.xyzOfVert, m.vertOfTria, thetasDegTria)
end

#test_compute_angle()

end
