#!/usr/bin/env julia
using SkyCoords
using BenchmarkTools

function myprintln(io, string)
    println(io, string)
    println(string)
end

io = open("julia_times.csv", "w")
println(io, "n,system,time")
for n in [1, 10, 100, 1000, 10_000, 100_000, 1_000_000]
    if n == 1
        c = ICRSCoords(2pi*rand(), pi*(rand() - 0.5))
        t1 = minimum(@benchmark(convert(GalCoords, $c)).times)/1e9
        myprintln(io, "$n,galactic,$t1")
        t2 = minimum(@benchmark(convert(FK5Coords{2000}, $c)).times)/1e9
        myprintln(io, "$n,fk5j2000,$t2")
        t3 = minimum(@benchmark(convert(FK5Coords{1975}, $c)).times)/1e9
        myprintln(io, "$n,fk5j1975,$t3")
    else
        c = [ICRSCoords(2pi*rand(), pi*(rand() - 0.5)) for i=1:n]
        t1 = minimum(@benchmark(convert(Vector{GalCoords{Float64}}, $c)).times)/1e9
        myprintln(io, "$n,galactic,$t1")
        t2 = minimum(@benchmark(convert(Vector{FK5Coords{2000, Float64}}, $c)).times)/1e9
        myprintln(io, "$n,fk5j2000,$t2")
        t3 = minimum(@benchmark(convert(Vector{FK5Coords{1975, Float64}}, $c)).times)/1e9
        myprintln(io, "$n,fk5j1975,$t3")
    end
    println()
end
close(io)
