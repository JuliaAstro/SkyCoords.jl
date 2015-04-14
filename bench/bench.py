#!/usr/bin/env python
from __future__ import division, print_function

import timeit

io = open("python_times.csv", "w")
nloops = 10
io.write("n,system,time\n")
for n in [1, 10, 100, 1000, 10000, 100000, 1000000]:
    # hack to determine number of loops
    nloops = 10 if n < 1000000 else 1
    setup = """
from numpy import pi
from numpy.random import rand
from astropy.coordinates import SkyCoord, FK5
n = {0:d}
if n == 1:
    c = SkyCoord(2.*pi*rand(n)[0], pi*rand(n)[0]-pi/2, unit=('rad', 'rad'))
else:
    c = SkyCoord(2.*pi*rand(n), pi*rand(n)-pi/2, unit=('rad', 'rad'))
""".format(n)
    t1 = timeit.repeat("c.galactic", setup, repeat=3, number=nloops)
    t2 = timeit.repeat("c.transform_to('fk5')", setup, repeat=3, number=nloops)
    t3 = timeit.repeat("c.transform_to(FK5(equinox='J1975'))", setup, repeat=3, number=nloops)
    line = """{n:d},galactic,{t1:f}
{n:d},fk5j2000,{t2:f}
{n:d},fk5j1975,{t3:f}
""".format(n=n, t1=min(t1)/nloops, t2=min(t2)/nloops, t3=min(t3)/nloops)
    io.write(line)
    print(line,)
io.close()
