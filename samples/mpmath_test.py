import mpmath
import mpfft
import pylab

mpmath.mp.dps=100
mpmath.mp.pretty=False

x = mpmath.linspace(0,mpmath.pi,100)
y = [mpmath.sin(xx) for xx in x ]
print(x[-1])
print(mpmath.pi)
print(y[-1])

a = mpmath.pi/mpmath.mpf(2)
b = mpmath.fdiv(mpmath.pi, mpmath.mpf(2))
print(mpmath.sin(a))
print(mpmath.sin(b))

x = [mpmath.mpc("0","0") for i in range(100)]
x[10] = mpmath.mpc("1","0")
z = mpfft.fft(x)
zreal = [z[i].real for i in range(len(z))]
zimag = [z[i].imag for i in range(len(z))]
#print z
print zreal
pylab.plot(zreal)
pylab.plot(zimag)
pylab.show()
