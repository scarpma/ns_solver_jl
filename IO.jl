module IO

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


function writeFieldsToFile(filename::String,x,y,u,v,p)
    vtkfile = vtk_grid(filename, x, y)
    vtkfile["U"] = vcat(u,v)
    vtkfile["p"] = staggeredToNormal(p)
    outfiles = vtk_save(vtkfile)
end

function writeSingleFieldToFile(filename::String,x,y,p,fieldName::String)
    vtkfile = vtk_grid(filename, x, y)
    vtkfile[fieldName] = staggeredToNormal(p)
    outfiles = vtk_save(vtkfile)
end

end
