#!/usr/bin/env julia
using CSV
using StatsPlots

df1 = CSV.read("julia_times.csv")
df2 = CSV.read("python_times.csv")

# rename systems
df1[:system] = ["jl_" * s for s in df1[:system]]
df2[:system] = ["py_" * s for s in df2[:system]]
df = vcat(df1, df2)

# Plot
@df df plot(:n, :time,
    group = :system,
    markershape = :circle,
    markersize = 4,
    markerstrokealpha = 0,
    linewidth = 1,
    alpha = 0.8,
    xlabel = "Number of coordinates",
    ylabel = "Time (s)",
    xscale = :log10,
    yscale = :log10,
    legend = :bottomright,
)
savefig("bench.png")
