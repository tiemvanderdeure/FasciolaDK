import AbstractGPs.KernelFunctions

struct DistanceMatrix{T} <: AbstractMatrix{T}
    A::Matrix{T}
end
Base.parent(A::DistanceMatrix) = A.A
Base.size(A::DistanceMatrix) = Base.size(parent(A))
Base.getindex(A::DistanceMatrix, i::Int) = getindex(parent(A), i)
Base.getindex(A::DistanceMatrix, I::Vararg{Int, N}) where N = getindex(parent(A), I...)
Base.setindex!(A::DistanceMatrix, v, I...) = setindex!(parent(A), v, I...)
Base.IndexStyle(::DistanceMatrix) = Base.IndexStyle(Matrix)
Base.similar(DA::DistanceMatrix) = DistanceMatrix(similar(parent(DA)))
function KernelFunctions.kernelmatrix(κ::KernelFunctions.SimpleKernel, A::DistanceMatrix)
    return map(x -> KernelFunctions.kappa(κ, x), A)
end
function KernelFunctions.kernelmatrix(κ::KernelFunctions.TransformedKernel, x::DistanceMatrix)
    return KernelFunctions.kernelmatrix(κ.kernel, DistanceMatrix(κ.transform.(x)))
end