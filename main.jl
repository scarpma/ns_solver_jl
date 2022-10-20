include("./IO.jl")
include("./FiniteDifference.jl")
using .IO
using .FiniteDifference

const Re = 100.
const dt = 1.0
const n1 = 50
const n2 = 50
const h1 = 1.0
const h2 = 1.0

x = collect(range(start=0., step=h2, length=n2))
y = collect(range(start=0., step=h1, length=n1))

xp = collect(range(start=-h2/2., step=h2, length=n2+1))
yp = collect(range(start=(n1-1)*h1+h1/2, step=-h1, length=n1+1))

u = zeros(Float32, (n1,n2))
v = zeros(Float32, (n1,n2))
p = zeros(Float32, (n1+1,n2+1))

CrankNich = FiniteDifference.CNO(n1*n2, dt, Re, h1, h2)
print(CrankNich[1,1])
# very inefficient. Gotta understand if the problem is in the
# implementation or it is normal. Is it possible that with
# a grid of 50x50 16Gb of RAM are not enough ? Very Strange !!!

IO.writeFieldsToFile("prova",x,y,u,v,p)
IO.writeSingleFieldToFile("prova_p",xp,yp,p,"p")
