module FiniteDifference

using BandedMatrices # BandedMatrix
using BlockBandedMatrices # Fill
using LazyArrays # Kron
using LinearAlgebra # I

function laplacian_operator(n, h1, h2)
    # 1d matrix for second order derivative
    D² = BandedMatrix(0 => Fill(-2,n), 1 => Fill(1,n-1), -1 => Fill(1,n-1))
    # 2d matrix for ∂_11
    D_11 = BandedBlockBandedMatrix(Kron(D²/h1^2., Eye(n)))
    # for 22
    D_22 = BandedBlockBandedMatrix(Kron(Eye(n), D²/h2^2.))
    # laplacian
    Δ = D_11 + D_22
    return Δ
end

function CNO(n, dt, Re, h1, h2)
    # Crank-Nicholson Operator in Matrix form
    Δ = @inline laplacian_operator(n, h1, h2)
    return I - (dt/(2*Re))*Δ
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

function test_LaplacianTriDiagSolve()
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

end
