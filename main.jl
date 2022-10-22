include("./IO.jl")
include("./FiniteDifference.jl")
using .IO
using .FiniteDifference

const T = Float64
const Re = 100.
const dt = 1.0
const n1 = 50
const n2 = 50
const h1 = 1.0
const h2 = 1.0


function main()

x = collect(range(start=0., step=h2, length=n2))
y = collect(range(start=0., step=h1, length=n1))

xp = collect(range(start=-h2/2., step=h2, length=n2+1))
yp = collect(range(start=(n1-1)*h1+h1/2, step=-h1, length=n1+1))

u = Field2D(T, (h1,h2), (n1,n2), 0.)
v = Field2D(T, (h1,h2), (n1,n2), 0.)
p = Field2D(T, (h1,h2), (n1+1,n2+1), 0.)

FiniteDifference.test_dd1()

#IO.writeFieldsToFile("prova",x,y,u,v,p)
#IO.writeSingleFieldToFile("prova_p",xp,yp,p,"p")

end


main()
