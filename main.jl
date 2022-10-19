include("./io.jl")
using .io

n1 = 50
n2 = 50
h1 = 1.0
h2 = 1.0

x = collect(range(start=0, stop=(n2-1)*h2, step=h1))
y = collect(range(start=(n1-1)*h1, stop=0, step=-h2))

xp = collect(range(start=-h1/2., stop=(n2)*h2, step=h1))
yp = collect(range(start=(n1)*h1+h2/2, stop=-h2/2, step=-h2))

u = zeros(Float32, (n1,n2))
v = zeros(Float32, (n1,n2))
p = zeros(Float32, (n1+1,n2+1))

function residual_u(u,v,p)
end
function residual_v(u,v,p)
end
function residual_p(u,v,p)
end

io.writeFieldsToFile(x,y,u,v,p)
