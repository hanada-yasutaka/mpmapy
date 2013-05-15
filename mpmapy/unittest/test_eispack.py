import mpmath
import HilbertSpace as H
mpmath.mp.dps=40
a = mpmath.mpc("0","1")
b = mpmath.mpc(mpmath.pi,"0")
c = mpmath.mpc(mpmath.pi,"0")
d = mpmath.mpc("0", "1")
z = mpmath.mpc("0","0")
#mat = mpmath.matrix([[a,b],[c,d]])
mat = mpmath.matrix([[a,b,z],[b,a,b],[z,b,a]])
Mat=H.Matrix(3,mat)
print Mat.data
Mat.eigen()
