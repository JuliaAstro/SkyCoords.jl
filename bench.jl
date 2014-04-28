using SkyCoords
using TimeIt

c1 = [ICRS(2pi*rand(), pi*(rand() - 0.5)) for i=1:100]

@timeit to_galactic(c1)
