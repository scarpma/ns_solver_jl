module IO
using Printf
using WriteVTK

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
