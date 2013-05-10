import HilbertSpace as H
import mpmath
import numpy
mpmath.mp.dps = 8
mpmath.mp.pretty = True
dim = 10
vec = H.Vector(dim)
print vec.data
xre = numpy.random.random(dim) -0.5
xim = numpy.random.random(dim) - 0.5
x = [mpmath.mpc(xre[i],xim[i]) for i in range(dim)]
y = xre + 1.j*xim
vec = H.Vector(dim, x)
print numpy.all(vec.data == x )
print vec.norm()
print vec.real()
