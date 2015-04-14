#!/usr/bin/env julia
using DataFrames
using Gadfly

df1 = readtable("julia_times.csv")
df2 = readtable("python_times.csv")

# rename systems
df1[:system] = ["jl_" * s for s in df1[:system]]
df2[:system] = ["py_" * s for s in df2[:system]]
df = vcat(df1, df2)

# Plot
p = plot(df, x="n", y="time", color="system", Geom.point, Geom.line,
         Scale.x_log10, Scale.y_log10, Guide.xlabel("# coordinates"),
         Guide.ylabel("time (s)"))
draw(PNG("bench.png", 7inch, 5.5inch), p)
