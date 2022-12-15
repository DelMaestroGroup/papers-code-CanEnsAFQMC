"""
    Canonical Ensemble Recursions

    # General Arguments
    Ns -> number of energy levels
    N -> desired particle number
    λ -> exponentiated spectrum, i.e., λ = exp(-βϵ)
"""

include("./quickmath.jl")

### Boson ###
function pf_recursion_boson(
    λ::AbstractVector{T}, N::Int64; 
    Ns = length(λ), 
    logZ = zeros(T, N + 1),
    n = zeros(T, Ns)
) where {T}
    """
    Borrmann recursion for bosons. Partition function along with occupation numbers are returned.
    """
    for i in 1 : N
        # rescale the term λ(1+n) to avoid numerical overflows
        logλn = log.(λ .* (1 .+ n))
        logλn_max = maximum(real(logλn))
        logλn .-= logλn_max
        λn = sum(exp.(logλn))

        # recursion is then expressed in the logarithmic form
        logZ[i + 1] = log(λn) + logλn_max - log(i) + logZ[i]
        n = exp.(logZ[i] - logZ[i + 1] + logλn_max .+ logλn)
    end

    # real part and phase are separated
    rlogZ = real(logZ[N + 1])
    sgnlogZ = exp(imag(logZ[N + 1])im)

    return rlogZ, sgnlogZ, n
end

### Fermion ###
function pf_recursion_fermion(
    λ::AbstractVector{T}, N::Int;
    isReal::Bool = false,
    Ns = length(λ),
    P::AbstractMatrix{Tp} = zeros(eltype(λ), N + 1, Ns)
) where {T, Tp}
    """
    Recursive calculation for the fermionic partition function
    """
    isReal || (λ = complex(λ))

    N == 0 && return convert(T, 0.0), 1.0
    if N == 1 
        logZ = log(sum(λ))
        rlogZ = real(logZ)
        sgnlogZ = exp(imag(logZ)im)

        return rlogZ, sgnlogZ
    elseif N == Ns
        logZ = sum(log.(λ))
        rlogZ = real(logZ)
        sgnlogZ = exp(imag(logZ)im)

        return rlogZ, sgnlogZ
    elseif N == Ns - 1
        logZ = sum(log.(λ)) - log(sum(λ))
        rlogZ = real(logZ)
        sgnlogZ = exp(imag(logZ)im)

        return rlogZ, sgnlogZ
    end
    
    # rescale spectrum
    expβμ = fermilevel(λ, N)
    expβϵμ = λ / expβμ
    
    # recursion in the probability form
    poissbino(expβϵμ, N, P=P)

    logP0 = -sum(log.(1 .+ expβϵμ))
    logZ = log(P[N+1, Ns]) - logP0 + N*log(expβμ)

    # real part and phase are separated
    rlogZ = real(logZ)
    sgnlogZ = exp(imag(logZ)im)

    return rlogZ, sgnlogZ
end

function pf_recursion(
    λ::Vector{T}, N::Int, ϵ::Float64;
    isReal::Bool = false,
    Ns = length(λ),
    P::AbstractMatrix{Tp} = zeros(eltype(λ), N + 1, Ns)
) where {T, Tp}
    """
    Recursive calculation of the partition function with low-temperature approximation

    ϵ -> tolerance for the partition function. Empirically, it can be directly used as
        the truncation threshold for the level occupancy, p.
    """
    isReal || (λ = complex(λ))

    expβμ = fermilevel(λ, N)
    expβϵμ = λ / expβμ

    p = expβϵμ ./ (1 .+ expβϵμ)
    Imp = 1 ./ (1 .+ expβϵμ)
    
    # determine the upper and lower bound of the active energy levels
    Nsu = 1
    while real(p[Nsu]) < ϵ
        Nsu += 1
    end
    Nsl = Ns
    while real(1 - p[Nsl]) < ϵ
        Nsl -= 1
    end

    # then the occupied levels and the active energy levels are picked up
    λocc = @view λ[Nsl + 1 : Ns]
    λf = @view expβϵμ[Nsu : Nsl]
    Nsf = length(λf)
    Nf = N - (Ns - Nsl)
    poissbino(λf, Nf, P=P, ν1 = (@view p[Nsu : Nsl]), ν2= (@view Imp[Nsu : Nsl]))

    logP0 = -sum(log.(1 .+ λf))
    logZ = log(P[Nf+1, Nsf]) - logP0 + Nf*log(expβμ) + sum(log.(λocc))

    rlogZ = real(logZ)
    sgnlogZ = exp(imag(logZ)im)

    return rlogZ, sgnlogZ
end

function occ_recursion(
    λ::Vector{T}, N::Int64;
    isReal::Bool = false,
    Ns::Int = length(λ),
    P::AbstractMatrix{Tp} = zeros(eltype(λ), Ns + 1, Ns)
) where {T, Tp}
    """
    Recursive calculation of the occupation number
    """
    isReal || (λ = complex(λ))

    N == 0 && return zeros(T, Ns)
    N == Ns && return ones(T, Ns)
    N == 1 && return λ / sum(λ)

    expβμ = fermilevel(λ, N)
    expβϵμ = λ / expβμ
    poissbino(expβϵμ, P=P)

    # num of energy levels below the Fermi level
    # use this formula to ensure complex conjugate pairs are in the same section
    N_below = sum(abs.(expβϵμ) .> 1)

    # separate recursions for occupancies of energy levels above/below the Fermi level
    @views n_above = occ_recursion_rescaled(Ns, N, expβϵμ[1 : Ns - N_below], P[:, Ns])
    @views n_below = occ_recursion_rescaled(Ns, N, expβϵμ[Ns - N_below + 1 : Ns], P[:, Ns], isReverse=true)
    # then concatenate
    n = vcat(n_above, n_below)

    return n
end

function occ_recursion_rescaled(
    Ns::Int64, N::Int64,
    expβϵμ::AbstractArray{Te}, P::AbstractArray{Tp};
    isReverse::Bool = false
) where {Te, Tp}
    """
    Level occupancy recursion using the rescaled spectrum
    """
    Ñs = length(expβϵμ)
    if !isReverse
        n = zeros(ComplexF64, N + 1, Ñs)
        @inbounds for i = 2 : N
            @views n[i + 1, :] = (P[i] / P[i + 1]) * expβϵμ .* (1 .- n[i, :])
            # Truncate the values that are smaller than 10^-10
            @views n[i + 1, :] = n[i + 1, :] .* (abs.(real(n[i + 1, :])) .> 1e-10)
        end

        return n[N + 1, :]
    else
        # num. of reverse recursions
        N_rev = Ns - N
        n = ones(ComplexF64, N_rev, Ñs)
        @inbounds for i = 1 : N_rev - 1
            @views n[i + 1, :] = (P[Ns - i + 1] / P[Ns - i]) * n[i, :] ./ expβϵμ
            # Truncate the values that are smaller than 10^-10
            @views n[i + 1, :] = n[i + 1, :] .* (abs.(real(n[i + 1, :])) .> 1e-10)
            @views n[i + 1, :] = 1 .- n[i + 1, :]
        end

        return n[N_rev, :]
    end
end
