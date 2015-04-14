#!/usr/bin/env julia
using SkyCoords
using TimeIt

io = open("julia_times.csv", "w")
println(io, "n,system,time")
for n in [1, 10, 100, 1000, 10_000, 100_000, 1_000_000]
    println("$n coordinates: ")
    if n == 1
        c = ICRSCoords(2pi*rand(), pi*(rand() - 0.5))
        t1 = @timeit convert(GalacticCoords, c)
        println(io, "$n,galactic,$t1")
        t2 = @timeit convert(FK5Coords{2000}, c)
        println(io, "$n,fk5j2000,$t2")
        t3 = @timeit convert(FK5Coords{1975}, c)
        println(io, "$n,fk5j1975,$t3")
    else
        c = [ICRSCoords(2pi*rand(), pi*(rand() - 0.5)) for i=1:n]
        t1 = @timeit convert(Vector{GalacticCoords}, c)
        println(io, "$n,galactic,$t1")
        t2 = @timeit convert(Vector{FK5Coords{2000}}, c)
        println(io, "$n,fk5j2000,$t2")
        t3 = @timeit convert(Vector{FK5Coords{1975}}, c)
        println(io, "$n,fk5j1975,$t3")
    end
end
close(io)
