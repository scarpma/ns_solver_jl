include("./IO.jl")
using .IO
using LinearAlgebra
using Printf

mutable struct SurfaceMesh{T}
    nv::Int64 # number of vertices
    ne::Int64 # number of edges
    nt::Int64 # number of triangles
    xyzOfVert::Array{T, 2}
    vertOfTria::Array{Int64, 2}
    vertOfEdge::Array{Int64, 2}
    edgeOfTria::Array{Int64, 2}
    lenOfEdges0::Array{T, 1}
    #vertOfVert::Array{Int64, 2}
    #nEdgeVert::Array{Int64, 1}
    #dirEdge::Array{T, 2}
    function SurfaceMesh(
        xyzOfVert::Union{Array{Float64, 2}, Array{Float32, 2}},
        vertOfTria::Array{Int64, 2},
        vertOfEdge::Array{Int64, 2},
        edgeOfTria::Array{Int64, 2},
    )
        nv = size(xyzOfVert, 2)
        ne = size(vertOfEdge, 2)
        nt = size(vertOfTria, 2)
        lenOfEdges0 = Array{eltype(xyzOfVert)}(undef, ne)
        m = new{eltype(xyzOfVert)}(nv, ne, nt, xyzOfVert, vertOfTria, vertOfEdge, edgeOfTria, lenOfEdges0)
        compute_distances!(m.lenOfEdges0, m)
        return m
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

function compute_cell_normals_and_areas!(normal, area, m)
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
    compute_cell_normals_and_areas!(normals, areas, m)
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

function expand_mesh(m, pnormals, A)
    m2 = deepcopy(m)
    for i=1:m.nv
        m2.xyzOfVert[:,i] = m.xyzOfVert[:,i] .+ A/0.027 * pnormals[:,i]
    end
    return m2
end

function compute_triaOfEdge(m)
    triaOfEdge = zeros(Int64, (2,m.ne))
    nTriaOfEdge = zeros(Int64, m.ne)
    for i=1:m.nt
        e1, e2, e3 = m.edgeOfTria[:,i]
        nTriaOfEdge[e1] += 1
        nTriaOfEdge[e2] += 1
        nTriaOfEdge[e3] += 1
        triaOfEdge[nTriaOfEdge[e1], e1] = i
        triaOfEdge[nTriaOfEdge[e2], e2] = i
        triaOfEdge[nTriaOfEdge[e3], e3] = i
    end
    return triaOfEdge, nTriaOfEdge
end

function compute_strains!(strains, m)
    d = zeros(eltype(strains), m.ne)
    compute_distances!(d, m)
    d0 = m.lenOfEdges0
    for i=1:m.ne
        strains[i] = (d[i] - d0[i]) / d[i]
    end
end

function internal_forces!(forces, m, areas, E2d)
    forces .= 0.0
    T = eltype(forces)
    potenergy::T = 0.0
    triaOfEdge, nTriaOfEdge = compute_triaOfEdge(m)
    d = zeros(eltype(forces), m.ne)
    compute_distances!(d, m)
    d0 = m.lenOfEdges0
    ke = zeros(Float64, m.ne)
    tmp = zeros(eltype(forces), 3)
    for i=1:m.ne
        t1, t2 = triaOfEdge[:,i]
        if (t1 != 0 && t2 != 0)
            ke[i] = E2d*(areas[t1]+areas[t2])/d[i]^2
        elseif (t1 != 0 && t2 == 0)
            ke[i] = E2d*(areas[t1])/d[i]^2
        else
            ke[i] = E2d*(areas[t2])/d[i]^2
        end
        v1, v2 = m.vertOfEdge[:,i]
        tmp .= -ke[i] * (d[i] - d0[i]) / d[i] .* (m.xyzOfVert[:,v1] .- m.xyzOfVert[:,v2])
        forces[:, v1] .+= +tmp[:]
        forces[:, v2] .+= -tmp[:]
        potenergy += ke[i] * (d[i] - d0[i])^2.
    end
    return potenergy
end

function velocity_verlet!(x,v,a,am,dt)
    x .= x + dt .* v + 0.5 .* am * dt.^2.
    v .= v + 0.5 .* (a + am) .* dt
end

compute_kenergy(v, mass) = 0.5 * mass * sum(v.^2.)

function compute_vertQuartetOfEdge(m)
    triaOfEdge, nTriaOfEdge = compute_triaOfEdge(m)
    vertQuartetOfEdge = zeros(Int64, (4, m.ne))
    for i=1:m.ne
        t1, t2 = triaOfEdge[:,i]
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
            vertQuartetOfEdge[:,i] .= v1, v2, v3, v4
        end
    end
    return vertQuartetOfEdge
end

function compute_angle!(theta, m::SurfaceMesh, normals)
    triaOfEdge, nTriaOfEdge = compute_triaOfEdge(m)
    for i=1:m.ne
        t1 = triaOfEdge[1,i]
        t2 = triaOfEdge[2,i]
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

function edgeArrayToTriaArray(m::SurfaceMesh, arr::Array{T,1})
    triaOfEdge, nTriaOfEdge = compute_triaOfEdge(m)
    T = eltype(m.xyzOfVert)
    newArr = zeros(T, m.nt)
    for i=1:m.ne
        t1, t2 = triaOfEdge[:,i]
        newArr[t1] += arr[i]
        if nTriaOfEdge[i] > 1
            newArr[t2] += arr[i]
        end
    end
    newArr = newArr ./ 3.
    return newArr
end

function test_compute_angle()
    T = Float64
    normals = zeros(T, (3,m.nt))
    areas = zeros(T, m.nt)
    thetas = zeros(T, m.ne)
    compute_cell_normals_and_areas!(normals, areas, m)
    compute_angle!(thetas, m, normals)
    thetasDegTria = edgeArrayToTriaArray(m, thetas) .* 360. / 2pi
    IO.writeSurface("provatheta", m.xyzOfVert, m.vertOfTria, thetasDegTria)
end

function main()
T = Float64
dt::T = 0.0001
tprint = 0.01
density = 1060*1.5
E = 1.e5 # young modulus [Pa]
thickness = 0.002 # thickness of the surface [m]
E2d = E*thickness
nsteps = 10000

coords, vertOfTria, vertOfEdge,
vertOfVert, edgeOfTria, nEdgeVert, 
midPointEdge, lenEdge, dirEdge = IO.read_gts_plus("./aorta1.gts", Float64)

m = SurfaceMesh(coords, vertOfTria, vertOfEdge, edgeOfTria)
triaOfVert, nTriaOfVert = compute_triaOfVert(m)
normals = zeros(T, (3,m.nt))
areas = zeros(T, m.nt)
compute_cell_normals_and_areas!(normals, areas, m)
@show sum(areas)
mass = density * sum(areas) * thickness / m.nv
@show mass
pnormals = cell_array_to_point_array(normals, m)
pnormals = pnormals ./ sqrt.(sum(pnormals.^2., dims=1))
IO.writeSurface("prova_orig", m.xyzOfVert, m.vertOfTria)


m2 = expand_mesh(m, pnormals, 0.0005)
forces = zeros(T, (3,m2.nv))
strains = zeros(T, m2.ne)
a = zeros(T,(3,m2.nv))
am = zeros(T,(3,m2.nv))
v = zeros(T,(3,m2.nv))
x = zeros(T,(3,m2.nv))
x = m2.xyzOfVert
potenergy = internal_forces!(forces, m2, areas, E2d)
kenergy = compute_kenergy(v, mass)

IO.writeSurface(@sprintf("prova_expanded_%5.5d",0), m2.xyzOfVert, m2.vertOfTria, forces, strains)
for i=1:nsteps
    @printf "---- t = %6.3f  i = %5d\n" dt*i  i
    a = forces ./ mass
    velocity_verlet!(x,v,a,am,dt)
    compute_strains!(strains, m2)
    @printf "  strain mean %.7e\n" sum(strains) / size(strains,1)
    tot_force = sum(forces, dims=2)
    @printf "  tot force   %10.3e %10.3e %10.3e\n" tot_force[1] tot_force[2] tot_force[3]
    @printf "  potenergy   %10.3e\n" potenergy
    @printf "  kenergy     %10.3e\n" kenergy
    @printf "  tot energy  %10.3e\n" kenergy + potenergy
    rms_force = sqrt.(sum(abs2.(v), dims=2)) ./ size(v,2)
    #@printf "  rms force %10.3e %10.3e %10.3e\n" rms_force[1] rms_force[2] rms_force[3]
    a_mean = sum(a, dims=2) / size(a, 2)
    #@printf "  mean acc %10.3e %10.3e %10.3e\n" a_mean[1] a_mean[2] a_mean[3]
    #@printf "  max acc %10.3e" maximum(a)
    #@printf "  min acc %10.3e\n" minimum(a)
    if mod(i, div(tprint, dt)) == 0
        IO.writeSurface(@sprintf("prova_expanded_%5.5d",i), m2.xyzOfVert, m2.vertOfTria, forces, strains)
    end
    potenergy = internal_forces!(forces, m2, areas, E2d)
    kenergy = compute_kenergy(v, mass)
    am = a
    m2.xyzOfVert = x
    println()
end
end

main()

