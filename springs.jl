include("./IO.jl")
include("./myMesh.jl")
using .IO
using .myMesh
using LinearAlgebra
using Printf

function expand_mesh(m, pnormals, A)
    m2 = deepcopy(m)
    for i=1:m.nv
        m2.xyzOfVert[:,i] = m.xyzOfVert[:,i] .+ A/0.027 * pnormals[:,i]
    end
    return m2
end

function compute_strains!(strains, m)
    d = zeros(eltype(strains), m.ne)
    compute_distances!(d, m)
    d0 = m.lenOfEdges0
    for i=1:m.ne
        strains[i] = (d[i] - d0[i]) / d[i]
    end
end

function internal_forces!(forces, m, E2d, ke, d, areas, thetas, normals)
    forces .= 0.0
    potenergy::eltype(forces) = 0.0
    compute_triaNormalAndArea!(normals, areas, m)
    compute_distances!(d, m)
    d0 = m.lenOfEdges0
    for i=1:m.ne
        t1, t2 = m.triaOfEdge[:,i]
        if (t1 != 0 && t2 != 0)
            ke[i] = E2d*(areas[t1]+areas[t2])/d[i]^2
        elseif (t1 != 0 && t2 == 0)
            ke[i] = E2d*(areas[t1])/d[i]^2
        else
            ke[i] = E2d*(areas[t2])/d[i]^2
        end
        v1, v2 = m.vertOfEdge[:,i]
        forces[:, v1] += -ke[i] * (d[i] - d0[i]) / d[i] .* (m.xyzOfVert[:,v1] .- m.xyzOfVert[:,v2])
        forces[:, v2] += -forces[:,v1]
        potenergy += ke[i] * (d[i] - d0[i])^2.
    end
    return potenergy
end

function velocity_verlet!(x,v,a,am,dt)
    x = x + dt .* v + 0.5 .* am .* dt^2.
    v = v + 0.5 .* (a + am) .* dt
end

compute_kenergy(v, mass) = 0.5 * mass * sum(v.^2.)

function main()
T = Float64
dt::T = 0.0001
tprint = 0.01
density = 1060*1.5
E = 1.e5 # young modulus [Pa]
thickness = 0.002 # thickness of the surface [m]
E2d = E*thickness
nsteps = 4

coords, vertOfTria, vertOfEdge,
vertOfVert, edgeOfTria, nEdgeVert, 
midPointEdge, lenEdge, dirEdge = IO.read_gts_plus("./aorta1.gts", Float64)

m = SurfaceMesh(coords, vertOfTria, vertOfEdge, edgeOfTria)
@show area = sum(m.triaAreas0)
mass = density * area * thickness / m.nv
@show mass
pnormals = myMesh.cell_array_to_point_array(m.triaNormals0, m)
pnormals = pnormals ./ sqrt.(sum(pnormals.^2., dims=1))
IO.writeSurface("prova_orig", m.xyzOfVert, m.vertOfTria)


m2 = expand_mesh(m, pnormals, 0.0005)
forces = zeros(T, (3,m2.nv))
ke = zeros(T, m.ne) 
strains = zeros(T, m2.ne)
a = zeros(T,(3,m2.nv))
am = zeros(T,(3,m2.nv))
v = zeros(T,(3,m2.nv))
x = zeros(T,(3,m2.nv))
d = zeros(eltype(forces), m.ne)
x = m2.xyzOfVert
areas = zeros(T,m2.nt)
normals = zeros(T,(3,m2.nt))
thetas = zeros(T,m2.ne)

for i=0:nsteps
    
    compute_triaNormalAndArea!(normals, areas, m2)
    compute_distances!(d, m2)
    #compute_angle!(thetas, m, normals)
    rms_v = sqrt.(sum(abs2.(v), dims=2)) ./ size(v,2)
    a_mean = sum(a, dims=2) / size(a, 2)
    potenergy = internal_forces!(forces, m2, E2d, ke, d, areas, thetas, normals)
    tot_force = sum(forces, dims=2)
    kenergy = compute_kenergy(v, mass)
    compute_strains!(strains, m2)
    a = forces ./ mass

    # write files if necessary
    if mod(i, div(tprint, dt)) == 0
        IO.writeSurface(@sprintf("prova_expanded_%5.5d",i), m2.xyzOfVert, m2.vertOfTria, forces, strains)
    end

    # compute stats
    @printf "---- t = %6.3f  i = %5d\n" dt*i  i
    @printf "  strain mean %.7e\n" sum(strains) / size(strains,1)
    @printf "  tot force   %10.3e %10.3e %10.3e\n" tot_force[1] tot_force[2] tot_force[3]
    @printf "  potenergy   %10.3e\n" potenergy
    @printf "  kenergy     %10.3e\n" kenergy
    @printf "  tot energy  %10.3e\n" kenergy + potenergy
    @printf "  rms velocity %10.3e %10.3e %10.3e\n" rms_v[1] rms_v[2] rms_v[3]
    #@printf "  mean acc %10.3e %10.3e %10.3e\n" a_mean[1] a_mean[2] a_mean[3]
    #@printf "  max acc %10.3e" maximum(a)
    #@printf "  min acc %10.3e\n" minimum(a)
    println()

    # integration
    velocity_verlet!(x,v,a,am,dt)
    m2.xyzOfVert = x
    am = a
end
end


main()
