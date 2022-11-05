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

function internal_forces!(forces, m, E2d, kb, ke, d, areas, thetas, normals, potenergy, tmp)
    forces[:,:] .= 0.0
    potenergy[1] = 0.0

    T = eltype(forces)
    b11::T = 0.
    b12::T = 1.
    b22::T = 0.
    beta::T = 0.
    a11 = zeros(T, 3)
    a12 = zeros(T, 3)
    a21 = zeros(T, 3)
    a22 = zeros(T, 3)
    a23 = zeros(T, 3)
    a31 = zeros(T, 3)
    a32 = zeros(T, 3)
    a33 = zeros(T, 3)
    a41 = zeros(T, 3)
    a42 = zeros(T, 3)

    # in plane springs
    compute_triaNormalAndArea!(normals, areas, m)
    compute_distances!(d, m)
    d0 = m.lenOfEdges0
    # update spring constants
    for i=1:m.ne
        t1, t2 = m.triaOfEdge[:,i]
        if (t1 != 0 && t2 != 0)
            ke[i] = E2d*(areas[t1]+areas[t2])/d[i]^2
        elseif (t1 != 0 && t2 == 0)
            ke[i] = E2d*(areas[t1])/d[i]^2
        else
            ke[i] = E2d*(areas[t2])/d[i]^2
        end
        # compute in plane spring forces
        v1, v2 = m.vertOfEdge[:,i]
        force = ( -ke[i] * (d[i] - d0[i]) / d[i] ) .* (m.xyzOfVert[:,v1] - m.xyzOfVert[:,v2])
        forces[:, v1] += force
        forces[:, v2] += -force
        potenergy[1] += ke[i] * (d[i] - d0[i])^2.
    end

    # bending springs
    thetas0 = m.angleOfEdges0
    compute_angle!(thetas, m, normals)
    for i=1:m.ne
        t1, t2 = m.triaOfEdge[:,i]
        # if edge is not on boundary
        if (t1 != 0 && t2 != 0)
            v1, v2, v3, v4 = m.vert1234OfEdge
            x1 = m.xyzOfVert[:,v1]
            x2 = m.xyzOfVert[:,v2]
            x3 = m.xyzOfVert[:,v3]
            x4 = m.xyzOfVert[:,v4]
            n1 = normals[:,t1]
            n2 = normals[:,t2]

            # compute coefficients b11 b12 b22 beta and cross products
            b11 = - cos(thetas[i])
            b22 = b11
            beta = kb * (sin(thetas[i])*cos(thetas0[i]) - cos(thetas[i])*sin(thetas0[i])) / sqrt(1. - cos(thetas[i])^2.)
            a11 = LinearAlgebra.cross(n1, x3 - x2)
            a12 = LinearAlgebra.cross(n2, x3 - x2)
            a21 = LinearAlgebra.cross(n1, x1 - x3)
            a22 = LinearAlgebra.cross(n1, x3 - x4) + LinearAlgebra.cross(n2, x1 - x3)
            a23 = LinearAlgebra.cross(n2, x3 - x4)
            a31 = LinearAlgebra.cross(n1, x2 - x1)
            a32 = LinearAlgebra.cross(n1, x3 - x2) + LinearAlgebra.cross(n2, x2 - x1)
            a33 = LinearAlgebra.cross(n2, x4 - x2)
            a41 = LinearAlgebra.cross(n1, x2 - x3)
            a42 = LinearAlgebra.cross(n2, x2 - x3)

            # compute bending spring forces on v1 v2 v3 v4
            forces[:, v1] += beta * (b11*a11 + b12*a12)
            forces[:, v2] += beta * (b11*a21 + b12*a22 + b22*a23)
            forces[:, v3] += beta * (b11*a31 + b12*a32 + b22*a33)
            forces[:, v4] += beta * (b11*a41 + b22*a42)
            potenergy[1] += kb * (1 - cos(thetas[i] - thetas0[i]))
        end
    end


end

function velocity_verlet!(
        x::Union{Array{Float64,2}, Array{Float32,2}},
        v::Union{Array{Float64,2}, Array{Float32,2}},
        a::Union{Array{Float64,2}, Array{Float32,2}},
        am::Union{Array{Float64,2}, Array{Float32,2}},
        dt::Union{Float64, Float32})
    x[:,:] = x + dt .* v + 0.5 .* am .* dt^2.
    v[:,:] = v + 0.5 .* (a + am) .* dt
end

function compute_kenergy!(kenergy, v, mass)
    kenergy[1] = 0.5 * mass * sum(v.^2.)
end

function main()
T = Float64
dt::T = 0.0001      # time step [s]
tprint = 0.0001       # save file time step [s]
density = 1060*1.5  # mass density of surface [kg / m^3]
E = 1.e5            # young modulus [Pa]
B = 1e-3             # bending modulus
thickness = 0.002   # thickness of the surface [m]
E2d = E*thickness
kb = B * 2. / sqrt(3.)
nsteps = 100

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
d = zeros(eltype(forces), m.ne)
areas = zeros(T,m2.nt)
normals = zeros(T,(3,m2.nt))
thetas = zeros(T,m2.ne)
potenergy = zeros(T, 1)
kenergy = zeros(T, 1)
tmp = zeros(T, 3)
tstrains = zeros(T, m2.nt)
dtheta = zeros(T, m2.nt)

for i=0:nsteps
    
    compute_triaNormalAndArea!(normals, areas, m2)
    compute_distances!(d, m2)
    #compute_angle!(thetas, m, normals)

    rms_v = sqrt.(sum(abs2.(v), dims=2) ./ size(v,2)) 
    a_mean = sum(a, dims=2) / size(a, 2)
    internal_forces!(forces, m2, E2d, kb, ke, d, areas, thetas, normals, potenergy, tmp)
    tot_force = sum(forces, dims=2)
    compute_kenergy!(kenergy, v, mass)
    compute_strains!(strains, m2)
    a = forces ./ mass

    # write files if necessary
    if mod(i, div(tprint, dt)) == 0
        myMesh.edgeArrayToTriaArray!(tstrains, m2, strains)
        myMesh.edgeArrayToTriaArray!(dtheta, m2, thetas - m.angleOfEdges0)
        IO.writeSurface(
            @sprintf("prova_expanded_%5.5d",i),
            m2.xyzOfVert, m2.vertOfTria,
            ["forces", "strains", "dtheta"],
            forces, tstrains, dtheta)
    end

    # compute stats
    @printf "---- t = %6.3f  i = %5d\n" dt*i  i
    @printf "  rms strain   %.7e\n" sqrt(sum(strains.^2.) / size(strains,1))
    @printf "  rms dtheta   %.7e\n" sqrt(sum(dtheta.^2.) / size(dtheta,1))
    @printf "  tot force    %10.3e %10.3e %10.3e\n" tot_force[1] tot_force[2] tot_force[3]
    @printf "  potenergy    %10.3e\n" potenergy[1]
    @printf "  kenergy      %10.3e\n" kenergy[1]
    @printf "  tot energy   %10.3e\n" kenergy[1] + potenergy[1]
    @printf "  rms velocity %10.3e %10.3e %10.3e\n" rms_v[1] rms_v[2] rms_v[3]
    #@printf "  mean acc %10.3e %10.3e %10.3e\n" a_mean[1] a_mean[2] a_mean[3]
    #@printf "  max acc %10.3e" maximum(a)
    #@printf "  min acc %10.3e\n" minimum(a)
    println()

    # integration
    velocity_verlet!(m2.xyzOfVert,v,a,am,dt)
    println("outside verlet ", m2.xyzOfVert[:,1])
    am = a
end
end


main()