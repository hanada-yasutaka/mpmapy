import mpMaps as Maps
import QmapSystem as QS
import mpmath
#import pylab
import time
import HilbertSpace as H
from QmapSystem import QmapSystem as Q
start=time.time()
dim=10

mpmath.mp.dps = 80
mpmath.mp.pretty=False
k=mpmath.mpf(1.0)
core = 1
verbose=True
map = Maps.StandardMap(k)
domain = [[0,1],[0,1]]
qmap = Q(map=map,dim=dim,type='U', domain=[[0,1],[0,1]])
evals,evecs=qmap.getEigen(core=core, verbose=verbose)
qmap.saveEigen()
qmap.saveHsm(core=core,verbose=verbose)
#print(evals)
#print(evecs[0])
#for vec in evecs:
#    vec.hsmrep(60,60)
#    exit()
    #print( isinstance(vec, QS.State))

exit()
"""
mat = H.Matrix(dim)
a = mpmath.mpc("0","1")
b = mpmath.mpc(mpmath.pi,"0")
mat.data =mpmath.matrix([[a,b],[b,a]])
mat.eigen()
#for eval in mat.evals:
#    print '###---- eigen values####'
print mpmath.pi
print 1/mpmath.sqrt(2)
print mat.evals.data
for evec in mat.evecs:
    print "###---- einge vectors ####"
    print evec.data
    
print mat.evecs[0].inner(mat.evecs[0])# == mpmath.mpc(1,0)
"""



U = QS.Unitary(map, dim)
U.setRange([[0,1],[0,1]])

#print U.scaleinfo.x
U.setMatrix(fft=True,verbose=True)
U.saveEigen(verbose=True)
exit()
eval, evecs = U.getEigen(verbose=True)

#pylab.plot(eval.real(), eval.imag(), 'o')
#pylab.show()

for i, vec in enumerate(evecs):
    pylab.plot(vec.abs2())
    pylab.semilogy()    
    pylab.savefig("eigen_vec_%d.png" % i)
    pylab.close()
    #pylab.show()
    #exit()

#mat =U.getMatrix()
#for i in range(3):
#    l = [mpmath.mpc("0","0") for i in range(dim)]
#    mat.data[i,:] = mpmath.matrix(l).T


#mat.eigen(varbose=True)

end = time.time()
print("main:", end - start, 'sec')