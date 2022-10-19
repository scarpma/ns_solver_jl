module io

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


function writeFieldsToFile(x,y,u,v,p)
    vtkfile = vtk_grid("my_vtk_file", x, y)
    vtkfile["U"] = vcat(u,v)
    vtkfile["p"] = staggeredToNormal(p)
    print(size(staggeredToNormal(p)))
    outfiles = vtk_save(vtkfile)
end

end
