using JLD, BenchmarkTools

include("./ce_recursion.jl")
include("./ce_projection.jl")

data = JLD.load("../data/mat_decomp_data_lowT.jld")
const Lx = collect(8:2:30)
const Pmat = zeros(ComplexF64, 900, 900)
const ϵ = 1e-8

t1_recursion = Float64[]
t2_recursion = Float64[]
t1_projection = Float64[]
t2_projection = Float64[]

for i in 1 : 14
    V = Lx[i] ^ 2
    N1 = div(V, 2)
    N2 = div(V, 10)
    λ = data["lambda"][i]

    # determine the minimum Fourier points
    Nft1 = 1
    logZ_exact = pf_projection(λ, N1)
    while abs(pf_projection(λ, N1, Nft = Nft1) - logZ_exact) > ϵ
        Nft1 += 1
    end

    Nft2 = 1
    logZ_exact = pf_projection(λ, N2)
    while abs(pf_projection(λ, N2, Nft = Nft2) - logZ_exact) > ϵ
        Nft2 += 1
    end

    tmp = @benchmark CanEnsAFQMC.pf_recursion($λ, $N1, $ϵ, P=$Pmat)
    push!(t1_recursion, mean(tmp.times) / 10^9)
    tmp = @benchmark CanEnsAFQMC.pf_recursion($λ, $N2, $ϵ, P=$Pmat)
    push!(t2_recursion, mean(tmp.times) / 10^9)

    tmp = @benchmark CanEnsAFQMC.pf_projection($λ, $N1, Nft=$Nft1)
    push!(t1_projection, mean(tmp.times) / 10^9)
    tmp = @benchmark CanEnsAFQMC.pf_projection($λ, $N2, Nft=$Nft2)
    push!(t2_projection, mean(tmp.times) / 10^9)
end

jldopen("../../data/benchmark_data_lowT.jld", "w") do file
    write(file, "t1_recursion", t1_recursion)
    write(file, "t2_recursion", t2_recursion)
    write(file, "t1_projection", t1_projection)
    write(file, "t2_projection", t2_projection)
end
