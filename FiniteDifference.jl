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

end
