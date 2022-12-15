using JLD, BenchmarkTools

include("./ce_recursion.jl")
include("./ce_projection.jl")

data = JLD.load("../../data/mat_decomp_data.jld")
const Lx = collect(8:2:30)
const Pmat = zeros(ComplexF64, 900, 900)

t1_recursion = Float64[]
t2_recursion = Float64[]
t1_projection = Float64[]
t2_projection = Float64[]

for i in 1 : 14
    V = Lx[i] ^ 2
    N1 = div(V, 2)
    N2 = div(V, 10)
    λ = data["lambda"][i]

    tmp = @benchmark pf_recursion_fermion($λ, $N1, P=$Pmat)
    push!(t1_recursion, mean(tmp.times) / 10^9)
    tmp = @benchmark pf_recursion_fermion($λ, $N2, P=$Pmat)
    push!(t2_recursion, mean(tmp.times) / 10^9)

    tmp = @benchmark pf_projection_fermion($λ, $N1)
    push!(t1_projection, mean(tmp.times) / 10^9)
    tmp = @benchmark pf_projection_fermion($λ, $N2)
    push!(t2_projection, mean(tmp.times) / 10^9)
end

jldopen("../data/benchmark_data.jld", "w") do file
    write(file, "t1_recursion", t1_recursion)
    write(file, "t2_recursion", t2_recursion)
    write(file, "t1_projection", t1_projection)
    write(file, "t2_projection", t2_projection)
end
