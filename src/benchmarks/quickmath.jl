"""
    Some math functions
"""
function fermilevel(expβϵ::Vector{T}, N::Int64) where {T<:Number}
    """
    Compute an approximate Fermi level
    """
    Ns = length(expβϵ)
    return sqrt(abs(expβϵ[Ns - N + 1] * expβϵ[Ns - N]))
end

function poissbino(
    ϵ::AbstractVector{T};
    Ns::Int64 = length(ϵ),
    ν1::AbstractVector{T} = ϵ ./ (1 .+ ϵ),
    ν2::AbstractVector{T} = 1 ./ (1 .+ ϵ),
    P::AbstractMatrix{Tp} = zeros(eltype(ϵ), Ns + 1, Ns)
) where {T<:Number, Tp<:Number}
    """
    A regularized version of the recursive calculation for the Poisson binomial
    distribution
    For details on how this is employed in the canonical ensemble calculations,
    see doi.org/10.1103/PhysRevResearch.2.043206

    # Argument
    ϵ -> eigenvalues
    """
    # Initialization
    P[1, 1] = ν2[1]
    P[2, 1] = ν1[1]
    # iteration over trials
    @inbounds for i = 2 : Ns
        P[1, i] = ν2[i] * P[1, i - 1]
        # iteration over number of successes
        for j = 2 : i + 1
            P[j, i] = ν2[i] * P[j, i - 1] + ν1[i] * P[j - 1, i - 1]
        end
    end

    return P
end

function poissbino(
    ϵ::AbstractVector{T}, N::Int64;
    Ns::Int64 = length(ϵ),
    ν1::AbstractVector{T} = ϵ ./ (1 .+ ϵ),
    ν2::AbstractVector{T} = 1 ./ (1 .+ ϵ),
    P::AbstractMatrix{Tp} = zeros(eltype(ϵ), N + 1, Ns)
) where {T<:Number, Tp<:Number}
    """
    Same Poisson binomial recursion scheme but stops at the Nth iteration
    """
    # Initialization
    P[1, 1] = ν2[1]
    P[2, 1] = ν1[1]
    # iteration over trials
    for i = 2 : Ns
        @inbounds P[1, i] = ν2[i] * P[1, i - 1]
        # iteration over number of successes
        for j = 2 : min(i + 1, N + 1)
            @inbounds P[j, i] = ν2[i] * P[j, i - 1] + ν1[i] * P[j - 1, i - 1]
        end
    end

    return P
end