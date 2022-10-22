module FiniteDifference

export Field2D

using BandedMatrices # BandedMatrix
using BlockBandedMatrices # Fill
using LazyArrays # Kron
using LinearAlgebra # I

using Plots

## function laplacian_operator(n, h1, h2)
##     # 1d matrix for second order derivative
##     D² = BandedMatrix(0 => Fill(-2,n), 1 => Fill(1,n-1), -1 => Fill(1,n-1))
##     # 2d matrix for ∂_11
##     D_11 = BandedBlockBandedMatrix(Kron(D²/h1^2., Eye(n)))
##     # for 22
##     D_22 = BandedBlockBandedMatrix(Kron(Eye(n), D²/h2^2.))
##     # laplacian
##     Δ = D_11 + D_22
##     return Δ
## end
## 
## function CNO(n, dt, Re, h1, h2)
##     # Crank-Nicholson Operator in Matrix form
##     Δ = @inline laplacian_operator(n, h1, h2)
##     return I - (dt/(2*Re))*Δ
## end

struct Field2D
    arr
    h1::AbstractFloat
    h2::AbstractFloat
    n1::Int64
    n2::Int64
    function Field2D(arr, (h1, h2))
        @assert ndims(arr) == 2
        n1 = size(arr,1)
        n2 = size(arr,2)
        new(arr, h1, h2, n1, n2)
    end
    function Field2D(T, (h1, h2), (n1,n2))
        arr = Array{T}(undef, (n1,n2))
        new(arr, h1, h2, n1, n2)
    end
    function Field2D(T, (h1, h2), (n1,n2), val)
        @assert typeof(val) == T
        arr = Array{T}(undef, (n1,n2))
        fill!(arr, val)
        new(arr, h1, h2, n1, n2)
    end
    Base.iterate(foo::Field2D, state = 1) = state > fieldcount(Field2D) ? nothing : (getfield(foo, state), state + 1)
end

function dd1!(u::Field2D, o::Field2D)

    arr, h1, h2, n1, n2 = u
    res, h1p, h2p, n1p, n2p = o
    @assert n1 == n1p
    @assert n2 == n2p
    @assert h1 == h1p
    @assert h2 == h2p

    #res = Array{eltype(arr)}(undef, (n1,n2))
    
    for j=1:n2
        res[1,j] = (arr[1,j] -2*arr[2,j] + arr[3,j]) / h1^2
        res[n1,j] = (arr[n1-2,j] -2*arr[n1-1,j] + arr[n1,j]) / h1^2
    end
    for i=2:n1-1
        for j=1:n2
            res[i,j] = (arr[i-1,j] -2*arr[i,j] + arr[i+1,j]) / h1^2
        end
    end
    return res
end

function test_dd1()
    f = Field2D(Float64, (0.1,0.1), (100,100))
    u, h1, h2, n1, n2 = f
    for j=1:n2
        for i=1:n1
            u[i,j] = cos(i*h1 + j*h2)
        end
    end
    p1 = heatmap(u, title="u")
    dd1u = dd1(f)
    p2 = heatmap(dd1u, title="d11 u")
    plot(p1, p2)
    gui()
end

function LaplacianTriDiagSolve!(v, d, h) # mutated argument in first position
    # particular case for toeplitz matrix with
    #   alpha = 1, beta = -2, gamma=1 : the 1d Laplacian discretization  
    # ai: lower diag    --> 1. / h
    # bi: diag          --> -2. / h
    # ci: upper diag    --> 1. / h
    # vi: variable (where the solution will be written)
    # di: RHS
    n = size(v,1)
    cp = Array{Float64}(undef, n-1)
    dp = Array{Float64}(undef, n)
    cp[1] = -1 / 2
    dp[1] = -d[1] / 2
    for i=2:n-1
        cp[i] = 1 / (-2 - cp[i-1])
        dp[i] = (d[i]-dp[i-1]) / (-2 - cp[i-1])
    end
    dp[n] = (d[n]-dp[n-1]) / (-2 - cp[n-1])
    #@show cp
    #@show dp
    v[n] = dp[n]*h^2
    for i=reverse(1:n-1)                            # backward substitution
        v[i] = dp[i]*h^2 - cp[i] * v[i+1]
    end
    return
end


function test_LaplacianTriDiagSolve1()
    n = 10000
    RHS = rand(Float64, (n))
    u = zeros(Float64, (n))
    FiniteDifference.LaplacianTriDiagSolve!(u, RHS, 1.0)
    @show u
    
    residual = 0.
    tmp = 0.

    tmp = (-2*u[1] + u[2])/h1 - RHS[1]
    println("Res[", 0, "]: ", tmp)
    residual += tmp
    for i=2:n-1
        tmp = (u[i-1]-2*u[i] + u[i+1])/h1 - RHS[i]
        println("Res[", i, "]: ", tmp)
        residual += tmp
    end
    tmp = (u[n-1] -2*u[n])/h1 - RHS[n]
    println("Res[", n, "]: ", tmp)
    residual += tmp
    @show residual
end
function test_LaplacianTriDiagSolve2()
    n = 1000
    f = Field2D(Float64, (1.0,1.0), (n,n), 0.0)
    RHS = rand(Float64, (n))
    u, h1 = f
    FiniteDifference.LaplacianTriDiagSolve!(@view(u[:,1]), RHS, 1.0)
    @show u
    
    residual = 0.
    tmp = 0.

    tmp = (-2*u[1] + u[2])/h1 - RHS[1]
    println("Res[", 0, "]: ", tmp)
    residual += tmp
    for i=2:n-1
        tmp = (u[i-1]-2*u[i] + u[i+1])/h1 - RHS[i]
        println("Res[", i, "]: ", tmp)
        residual += tmp
    end
    tmp = (u[n-1] -2*u[n])/h1 - RHS[n]
    println("Res[", n, "]: ", tmp)
    residual += tmp
    @show residual
end

end
