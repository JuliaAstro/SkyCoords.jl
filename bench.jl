#!/usr/bin/env julia
using SkyCoords
using TimeIt

for n in [10, 1000, 100000]
    print("$n coordinates: ")
    c = [ICRS(2pi*rand(), pi*(rand() - 0.5)) for i=1:n]
    @timeit to_galactic(c)
end
