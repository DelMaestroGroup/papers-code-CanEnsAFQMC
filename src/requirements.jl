using Pkg

dependencies = [
    "IJulia",
    "JLD",
    "Plots",
    "LaTeXStrings", 
    "DataFrames", 
    "GLM",
    "BenchmarkTools",
    "Measurements"
]

Pkg.add(dependencies)
