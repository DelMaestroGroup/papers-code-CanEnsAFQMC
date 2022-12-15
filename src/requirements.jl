using Pkg

dependencies = [
    "IJulia",
    "JLD",
    "Plots",
    "LaTeXStrings", 
    "DataFrames", 
    "GLM",
    "BenchmarkTools"
]

Pkg.add(dependencies)