"""
    Canonical Ensemble Projections

    # General Arguments
    Ns -> number of energy levels
    N -> desired particle number
    λ -> exponentiated spectrum, i.e., λ := exp(-βϵ)
    expiφ -> discrete Fourier frequencies
"""

include("./quickmath.jl")

### Boson ###
function pf_projection_boson(
    λ::AbstractVector{T}, N::Int64, M::Int64;
    # chemical potential is by default set to 0
    expβμ::Float64 = 1.0,
    iφ = im * [2 * π * m / M for m = 1 : M],
    expiφ = exp.(iφ)
) where T
    """
    Projection algorithm for bosonic partition functions
    """
    expiφβμ = expiφ / expβμ

    η = [sum(log.(1 .- expiφβμ[m] * λ)) for m = 1 : M]

    logZm = N*(-iφ .+ log(expβμ)) .- η
    logZm_max = maximum(real(logZm))
    Z = sum(exp.(logZm .- logZm_max))

    logZ = log(Z) + logZm_max - log(M)

    # real part and phase are separated
    rlogZ = real(logZ[N + 1])
    sgnlogZ = exp(imag(logZ[N + 1])im)

    return rlogZ, sgnlogZ
end

### Fermion ###
function pf_projection_fermion(
    λ::Vector{T}, N::Int; 
    returnFull::Bool = false,
    # number of Fourier points
    Nft = length(λ),
    iφ = im * [2 * π * m / Nft for m = 1 : Nft],
    expiφ = exp.(iφ)
) where T
    """
    Projection algorithm for fermionic partition functions
    """
    expβμ = fermilevel(λ, N)
    expiφmβμ = expiφ / expβμ

    η = [sum(log.(1 .+ expiφmβμ[m] * λ)) for m = 1 : Nft]

    logZm = N*(-iφ .+ log(expβμ)) .+ η

    # rescale all components by dividing the maximum value
    logZm_max = maximum(real(logZm))
    Z = sum(exp.(logZm .- logZm_max))

    logZ = log(Z) + logZm_max - log(Nft)

    # real part and phase are separated
    rlogZ = real(logZ[N + 1])
    sgnlogZ = exp(imag(logZ[N + 1])im)

    returnFull || return rlogZ, sgnlogZ
    return expβμ, expiφmβμ, η, logZ
end

function occ_projection_fermion(
    λ::Vector{T}, N::Int;
    Ns = length(λ),
    # number of Fourier points
    Nft = length(λ),
    iφ = im * [2 * π * m / Nft for m = 1 : Nft],
    expiφ = exp.(iφ),
    n = zeros(ComplexF64, Ns)
) where T
    """
    Projection algorithm for occupation numbers
    """
    expβμ, expiφmβμ, η, logZ = pf_projection_fermion(λ, N, returnFull=true, Nft=Nft, iφ=iφ, expiφ=expiφ)

    @inbounds for i in 1 : Ns
        ñk = expiφmβμ * λ[i] ./ (1 .+ expiφmβμ * λ[i])
        logñ = -iφ * N .+ log.(ñk) .+ η
        logn = logñ .+ log(expβμ) * N .- logZ
        n[i] = sum(exp.(logn)) / Nft
    end

    return n
end
