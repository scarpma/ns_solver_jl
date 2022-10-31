module FiniteDifference

global const DEBUG = true

export Field2D, dassert

using Printf
using Plots
include("./IO.jl")
using .IO

## using BandedMatrices # BandedMatrix
## using BlockBandedMatrices # Fill
## using LazyArrays # Kron
## using LinearAlgebra # I
## 
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
    h1::Float64
    h2::Float64
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

macro dassert(exp)
    if DEBUG
        return esc(:(@assert $exp))
    end
end

function d1!(res, arr, h)
    n1 = size(arr, 1)
    n2 = size(arr, 2)
    twice_h = 2*h
    for j=1:n2
        # x1 left border
        #res[1,j] = (res[2,j] - res[1,j]) / h1 # forward differentiation
        res[1,j] = (1.5*arr[1,j] - 2*arr[2,j] + 0.5*arr[3,j]) / h # forward differentiation
        # x1 right border
        #res[n1,j] = (res[n1,j] - res[n1-1,j]) / h1 # backward differentiation
        res[n1,j] = (1.5*arr[n1,j] - 2*arr[n1-1,j] + 0.5*arr[n1-2,j]) / h # backward differentiation
    end
    for j=1:n2
        for i=2:n1-1
            # inner part of grid (no borders)
            res[i,j] = (arr[i+1,j] - arr[i-1,j]) / twice_h # central differentiation
        end
    end
end
function d1!(res_field::Field2D, arr_field::Field2D)
    arr, h1, h2, n1, n2 = arr_field
    res, = res_field
    d1!(res, arr, h1)
end
function d1(arr_field::Field2D)
    arr, h1, h2, n1, n2 = arr_field
    res = Array{eltype(arr)}(undef, (n1,n2))
    d1!(res, arr, h1)
    return Field2D(res, (h1,h2))
end
function d2!(res, arr, h)
    n1 = size(arr, 1)
    n2 = size(arr, 2)
    twice_h = 2*h
    for i=1:n1
        # x1 left border
        #res[1,j] = (res[2,j] - res[1,j]) / h1 # forward differentiation
        res[i,1] = (1.5*arr[i,1] - 2*arr[i,2] + 0.5*arr[i,3]) / h # forward differentiation
        # x1 right border
        #res[n1,j] = (res[n1,j] - res[n1-1,j]) / h1 # backward differentiation
        res[i,n2] = (1.5*arr[i,n2] - 2*arr[i,n2-1] + 0.5*arr[i,n2-2]) / h # backward differentiation
    end
    for j=2:n2-1
        for i=1:n1
            # inner part of grid (no borders)
            res[i,j] = (arr[i,j+1] - arr[i,j-1]) / twice_h # central differentiation
        end
    end
end
function d1!(res_field::Field2D, arr_field::Field2D)
    arr, h1, h2, n1, n2 = arr_field
    res, = res_field
    d1!(res, arr, h1)
end
function d2!(res_field::Field2D, arr_field::Field2D)
    arr, h1, h2, n1, n2 = arr_field
    res, = res_field
    d2!(res, arr, h1)
end
function d1(arr_field::Field2D)
    arr, h1, h2, n1, n2 = arr_field
    res = Array{eltype(arr)}(undef, (n1,n2))
    d1!(res, arr, h1)
    return Field2D(res, (h1,h2))
end
function d2(arr_field::Field2D)
    arr, h1, h2, n1, n2 = arr_field
    res = Array{eltype(arr)}(undef, (n1,n2))
    d2!(res, arr, h1)
    return Field2D(res, (h1,h2))
end

function nonLinear1!(H1_field, u1_field, u2_field, d1u1_field, d2u1_field)
    # compute (non-linear term)_1 = H1 = (u1*d1+u2*d2)u1 = u1*d1u1 + u2*d2u1
    # compute d1u1 and d2u1
    d1!(d1u1_field, u1_field)
    d2!(d2u1_field, u1_field)

    u1, h1, h2, n1, n2 = u1_field
    u2, = u2_field
    H1, = H1_field
    d1u1, = d1u1_field
    d2u1, = d2u1_field

    for j=1:n2
        for i=1:n1
            H1[i,j] = u1[i,j] * d1u1[i,j] + u2[i,j] * d2u1[i,j]
        end
    end
end

function nonLinear2!(H2_field, u1_field, u2_field, d1u2_field, d2u2_field)
    # compute (non-linear term)_2 = (u1*d1+u2*d2)u2 = u1*d1u2 + u2*d2u2
    # compute d2u2 and d2u2
    d1!(d1u2_field, u2_field)
    d2!(d2u2_field, u2_field)

    u1, h1, h2, n1, n2 = u1_field
    u2, = u2_field
    H2, = H2_field
    d1u2, = d1u2_field
    d2u2, = d2u2_field

    for j=1:n2
        for i=1:n1
            H2[i,j] = u1[i,j] * d1u2[i,j] + u2[i,j] * d2u2[i,j]
        end
    end

end

function dd1!(res, arr, h)
    n1 = size(arr, 1)
    n2 = size(arr, 2)
    for j=1:n2
        # x1 left border
        #res[1,j] = (arr[1,j] -2*arr[2,j] +arr[3,j]) / h^2 # forward differentiation
        res[1,j] = (2*arr[1,j] -5*arr[2,j] + 4*arr[3,j] -arr[4,j]) / h^2 # forward differentiation
        # x1 right border
        #res[n1,j] = (arr[n1-2,j] -2*arr[n1-1,j] + arr[n1,j]) / h^2 # backward differentiation
        res[n1,j] = (-arr[n1-3,j] +4*arr[n1-2,j] -5*arr[n1-1,j] + 2*arr[n1,j]) / h^2 # backward differentiation
    end
    for i=2:n1-1
        for j=1:n2
            # inner part of grid (no borders)
            res[i,j] = (arr[i-1,j] -2*arr[i,j] + arr[i+1,j]) / h^2 # central differentiation
        end
    end
    return
end
function dd2!(res, arr, h)
    n1 = size(arr, 1)
    n2 = size(arr, 2)
    for i=1:n1
        # x2 left border
        #res[i,1] = (arr[i,1] -2*arr[i,2] + arr[i,3]) / h^2 # forward differentiation
        res[i,1] = (2*arr[i,1] -5*arr[i,2] +4*arr[i,3] -arr[i,4]) / h^2 # forward differentiation
        # x2 right border
        #res[i,n2] = (arr[i,n2-2] -2*arr[i,n2-1] + arr[i,n2]) / h^2 # backward differentiation
        res[i,n2] = (-arr[i,n2-3] +4*arr[i,n2-2] -5*arr[i,n2-1] + 2*arr[i,n2]) / h^2 # backward differentiation
    end
    for j=2:n2-1
        for i=1:n1
            # inner part of grid (no borders)
            res[i,j] = (arr[i,j-1] -2*arr[i,j] + arr[i,j+1]) / h^2 # central differentiation
        end
    end
    return
end
function dd1!(o::Field2D, u::Field2D)
    arr, h1, h2, n1, n2 = u
    res, h1p, h2p, n1p, n2p = o
    @dassert(n1 == n1p)
    @dassert(n2 == n2p)
    @dassert(h1 == h1p)
    @dassert(h2 == h2p)
    dd1!(res, arr, h1)
    return
end
function dd1!(u::Field2D)
    arr, h1, h2, n1, n2 = u
    res = Array{eltype(arr)}(undef, (n1,n2))
    dd1!(res, arr, h1)
    return Field2D(res, (h1,h2))
end
function dd2!(o::Field2D, u::Field2D)
    arr, h1, h2, n1, n2 = u
    res, h1p, h2p, n1p, n2p = o
    @dassert(n1 == n1p)
    @dassert(n2 == n2p)
    @dassert(h1 == h1p)
    @dassert(h2 == h2p)
    dd2!(res, arr, h2)
    return
end
function dd2!(u::Field2D)
    arr, h1, h2, n1, n2 = u
    res = Array{eltype(arr)}(undef, (n1,n2))
    dd2!(res, arr, h2)
    return Field2D(res, (h1,h2))
end


## TESTS
function compute_residual(sol, an_sol)
    n1 = size(sol, 1)
    n2 = size(sol, 2)
    tmp = 0.
    res = 0.
    max_res = 0.
    for j=1:n2
        for i=1:n1
            tmp = abs(sol[i,j] - an_sol[i,j])
            res += tmp
            if (tmp > max_res)
                max_res = tmp 
            end
        end
    end
    @printf "%10s %.3e\n" "sum res:" res
    res = res / (n1*n2)
    @printf "%10s %.3e\n" "mean res:" res
    @printf "%10s %.3e\n" "max res:" max_res
    return res
end
function getTestFunc1(n1,n2)
    res = Field2D(Float64, (2*pi/n1,2*pi/n2), (n1,n2))
    i, h1, h2, n1, n2 = res
    for jj=1:n2
        for ii=1:n1
            i[ii,jj] = cos(ii*h1 + jj*h2)
        end
    end
    return res
end
function getTestFunc2(n1,n2)
    res = Field2D(Float64, (2*pi/n1,2*pi/n2), (n1,n2))
    i, h1, h2, n1, n2 = res
    for jj=1:n2
        for ii=1:n1
            #i[ii,jj] = cos(ii*h1) * sin(jj*h2)
            i[ii,jj] = ((ii-1)*h1)^3 + ((jj-1)*h2)^3
        end
    end
    x1 = Array{Float64}(undef, (n1))
    x2 = Array{Float64}(undef, (n2))
    for ii=1:n1
        x1[ii] = (ii-1)*h1
    end
    for ii=1:n2
        x2[ii] = (ii-1)*h2
    end
    return res, x1, x2
end
function test_d1(n)
    n1 = n
    n2 = n
    input, x1, x2 = getTestFunc2(n1,n2)
    i, h1, h2 = input
    output = Field2D(Float64, (h1,h2), (n1, n2))
    o, = output

    #p1 = heatmap(i, title="input", aspect_ratio=:equal, xlabel="x_1", ylabel="x_2")
    d1!(output, input)
    IO.writeFieldsToFile("prova_d1", x1, x2, input.arr, output.arr)
    res = compute_residual(output.arr, [3*(i-1)^2*h1^2 for i in 1:n1, j in 1:n2])
    #p2 = heatmap(o, title="d11 input", aspect_ratio=:equal, xlabel="x_1", ylabel="x_2")
    #p = plot(p1, p2)
    #display(p)
    return res
end
function test_d2(n)
    n1 = n
    n2 = n
    input, x1, x2 = getTestFunc2(n1,n2)
    i, h1, h2 = input
    output = Field2D(Float64, (h1,h2), (n1, n2))
    o, = output

    #p1 = heatmap(i, title="input", aspect_ratio=:equal, xlabel="x_1", ylabel="x_2")
    d2!(output, input)
    IO.writeFieldsToFile("prova_d2", x1, x2, input.arr, output.arr)
    res = compute_residual(output.arr, [3*(j-1)^2*h1^2 for i in 1:n1, j in 1:n2])
    #p2 = heatmap(o, title="d11 input", aspect_ratio=:equal, xlabel="x_1", ylabel="x_2")
    #p = plot(p1, p2)
    #display(p)
    return res
end
function test_dd1(n)
    n1 = n
    n2 = n
    input, x1, x2 = getTestFunc2(n1,n2)
    i, h1, h2 = input
    output = Field2D(Float64, (h1,h2), (n1, n2))
    o, = output
    res = 0.

    #p1 = heatmap(i, title="input", aspect_ratio=:equal, xlabel="x_1", ylabel="x_2")
    dd1!(output, input)
    IO.writeFieldsToFile("prova_dd1", x1, x2, input.arr, output.arr)
    res = compute_residual(output.arr, [6*(i-1)*h1 for i in 1:n1, j in 1:n2])

    #p2 = heatmap(o, title="d11 input", aspect_ratio=:equal, xlabel="x_1", ylabel="x_2")
    #p = plot(p1, p2)
    #display(p)
    return res
end
function test_dd2(n)
    n1 = n
    n2 = n
    input, x1, x2 = getTestFunc2(n1,n2)
    i, h1, h2 = input
    output = Field2D(Float64, (h1,h2), (n1, n2))
    o, = output
    res = 0.

    #p1 = heatmap(i, title="input", aspect_ratio=:equal, xlabel="x_1", ylabel="x_2")
    dd2!(output, input)
    IO.writeFieldsToFile("prova_dd2", x1, x2, input.arr, output.arr)
    res = compute_residual(output.arr, [6*(j-1)*h2 for i in 1:n1, j in 1:n2])

    #p2 = heatmap(o, title="d11 input", aspect_ratio=:equal, xlabel="x_1", ylabel="x_2")
    #p = plot(p1, p2)
    #display(p)
    return res
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
