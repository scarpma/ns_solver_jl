include("./IO.jl")
include("./myMesh.jl")
using .IO
using .myMesh
using LinearAlgebra: dot
using Printf

function expand_mesh(m, pnormals, A)
    m2 = deepcopy(m)
    for i=1:m.nv
        m2.xyzOfVert[:,i] = m.xyzOfVert[:,i] .+ A/0.027 * pnormals[:,i]
    end
    return m2
end

function cross!(res, a, b)
    res[1] = a[2] * b[3] - a[3] * b[2]
    res[2] = a[3] * b[1] - a[1] * b[3]
    res[3] = a[1] * b[2] - a[2] * b[1]
end

function internal_forces!(forces, m, E2d, kb, visc, ke, d, areas, thetas, normals, v, potenergy, tmp)
    T = eltype(forces)
    forces[:,:] .= 0.0
    potenergy[1] = 0.0
    
    # for bending forces computations
    b11::T   = 0
    b12::T   = 0
    b22::T   = 0
    beta::T  = 0
    modn1::T = 0
    modn2::T = 0
    a21      = Vector{T}(undef, 3)
    a31      = Vector{T}(undef, 3)
    a32      = Vector{T}(undef, 3) 
    a13      = Vector{T}(undef, 3) 
    a34      = Vector{T}(undef, 3)
    a42      = Vector{T}(undef, 3)
    a23      = Vector{T}(undef, 3)
    a24      = Vector{T}(undef, 3)
    n1       = Vector{T}(undef, 3)
    n2       = Vector{T}(undef, 3)
    tvec1    = Vector{T}(undef, 3)
    tvec2    = Vector{T}(undef, 3)
    tvec3    = Vector{T}(undef, 3)
    tvec4    = Vector{T}(undef, 3)
    vel1     = Vector{T}(undef, 3)
    vel2     = Vector{T}(undef, 3)

    d0 = m.lenOfEdges0
    thetas0 = m.angleOfEdges0
    compute_angle!(thetas, m, normals)

    for i=1:m.ne
        t1, t2 = m.triaOfEdge[:,i]

        # in plane springs
        # update spring constants
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

        # bending springs
        # if edge is not on boundary
        if (t1 != 0 && t2 != 0)
            v1, v2, v3, v4 = m.vert1234OfEdge[:,i]
            x1 = m.xyzOfVert[:,v1]
            x2 = m.xyzOfVert[:,v2]
            x3 = m.xyzOfVert[:,v3]
            x4 = m.xyzOfVert[:,v4]
            n1 = normals[:,t1]
            n2 = normals[:,t2]

            a21[:] = x2 - x1
            a31[:] = x3 - x1
            a32[:] = x3 - x2
            a13[:] = x1 - x3
            a34[:] = x3 - x4
            a42[:] = x4 - x2
            a23[:] = x2 - x3
            a24[:] = x2 - x4

            cross!(n1, a21, a31)
            cross!(n2, a34, a24)
            modn1 = sqrt(sum(n1.^2.))
            modn2 = sqrt(sum(n2.^2.))

            # compute coefficients b11 b12 b22 beta and cross products
            beta = kb * (sin(thetas[i])*cos(thetas0[i]) - (dot(n1,n2)/(modn1*modn2))*sin(thetas0[i])) /
                #sqrt(1. - (dot(n1,n2)/(modn1*modn2))^2.)
                sqrt(1. - cos(thetas[i])^2.)
            b11 = -beta * (dot(n1,n2)/(modn1*modn2))/(modn1^2.)
            b12 = +beta / (modn1 * modn2)
            b22 = -beta * (dot(n1,n2)/(modn1*modn2))/(modn2^2.)

            # v1
            cross!(tvec1, n1, a32)
            cross!(tvec2, n2, a32)
            forces[:,v1] = forces[:,v1] + b11.*tvec1 + b12.*tvec2

            # v2
            cross!(tvec1, n1, a13)
            cross!(tvec2, n1, a34)
            cross!(tvec3, n2, a13)
            cross!(tvec4, n2, a34)
            forces[:,v2] = forces[:,v2] + b11.*tvec1 + b12.*(tvec2+tvec3) + b22.*tvec4

            # v3
            cross!(tvec1, n1, a21)
            cross!(tvec2, n1, a42)
            cross!(tvec3, n2, a21)
            cross!(tvec4, n2, a42)
            forces[:,v3] = forces[:,v3] + b11.*tvec1 + b12.*(tvec2+tvec3) + b22.*tvec4

            # v4
            cross!(tvec1, n1, a23)
            cross!(tvec2, n2, a23)
            forces[:,v4] = forces[:,v4] + b12.*tvec1 + b22.*tvec2

            potenergy[1] += kb * (1 - cos(thetas[i] - thetas0[i]))

            # damping
            v1, v2 = m.vertOfEdge[:,i]
            vel1[:] = v[:, v1]
            vel2[:] = v[:, v2]
            forces[:, v1] = forces[:, v1] - visc .* (vel1 - vel2)
            
            if (any(isnan, forces))
                println("IS NAN, i=", i)
                exit()
            end
 
        end
    end


end

function velocity_verlet!(
        x::Union{Array{Float64,2}, Array{Float32,2}},
        v::Union{Array{Float64,2}, Array{Float32,2}},
        a::Union{Array{Float64,2}, Array{Float32,2}},
        am::Union{Array{Float64,2}, Array{Float32,2}},
        dt::Union{Float64, Float32})
    for i=1:size(x,2)
        x[:,i] = x[:,i] + dt .* v[:,i] + 0.5 .* am[:,i] .* dt^2.
        v[:,i] = v[:,i] + 0.5 .* (a[:,i] + am[:,i]) .* dt
    end
end

function compute_kenergy!(kenergy, v, mass)
    kenergy[1] = 0.5 * mass * sum(v.^2.)
end

function main()
T = Float64
dt::T = 2e-5        # time step [s]
tprint = 1e-3       # save file time step [s]
density = 1.5       # mass density of surface [kg / m^3]
E = 719             # young modulus [Pa]
#E = 0.              # young modulus [Pa]
B = 0.325           # bending modulus
visc = 0.35         # damping coefficient
thickness = 0.0741  # thickness of the surface [m]
E2d = E*thickness
kb = B * 2. / sqrt(3.)
nsteps = 10000

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

m2 = expand_mesh(m, pnormals, 0.0015)
forces = zeros(T, (3,m2.nv))
ke = zeros(T, m.ne) 
strains = zeros(T, m2.ne)
a = zeros(T,(3,m2.nv))
am = zeros(T,(3,m2.nv))
v = zeros(T,(3,m2.nv))
d = zeros(eltype(forces), m.ne)
d0 = deepcopy(m.lenOfEdges0)
areas = zeros(T,m2.nt)
normals = zeros(T,(3,m2.nt))
thetas = zeros(T,m2.ne)
potenergy = zeros(T, 1)
kenergy = zeros(T, 1)
tmp = zeros(T, 3)
tstrains = zeros(T, m2.nt)
dtheta = zeros(T, m2.nt)

for i=0:nsteps
    
    compute_distances!(d, m2)
    compute_triaNormalAndArea!(normals, areas, m2)
    compute_angle!(thetas, m2, normals)

    rms_v = sqrt.(sum(abs2.(v), dims=2) ./ size(v,2)) 
    a_mean = sum(a, dims=2) / size(a, 2)
    internal_forces!(forces, m2, E2d, kb, visc, ke, d, areas, thetas, normals, v, potenergy, tmp)
    tot_force = sum(forces, dims=2)
    compute_kenergy!(kenergy, v, mass)
    strains[:] = (d[:] - d0[:]) ./ d[:]
    a[:,:] = forces ./ mass

    # write files if necessary
    if mod(i*dt, tprint) < dt
        myMesh.edgeArrayToTriaArray!(tstrains, m2, strains)
        myMesh.edgeArrayToTriaArray!(dtheta, m2, thetas - m.angleOfEdges0)
        IO.writeSurface(
            @sprintf("prova_expanded_%5.5d",Int64(round(i*dt/tprint))),
            m2.xyzOfVert, m2.vertOfTria,
            ["forces", "strains", "dtheta"],
            forces, tstrains, dtheta./(360/2pi))
    end

    # compute stats
    @printf "---- t = %7.4f  i = %5d\n" dt*i  i
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
    am[:,:] = a[:,:]
end
end


main()
