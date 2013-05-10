import unittest
import random
import mpmath
import numpy
import HilbertSpace as H


class TestQmap(unittest.TestCase):
    def setUp(self):
        mpmath.mp.dps=40
        mpmath.mp.pretty=True	


    def _testSetRange(self):
        qmin = random.randint(-10,0)
        pmin = random.randint(-10,0)

        qmax= random.randint(1,10)
        pmax = random.randint(1,10)
        self.S.setRange([[qmin,qmax],[pmin,pmax]])
        print(self.S.getRange())
        print(self.S.getX())
        print(self.S.getVolume())
        print(self.S.getPlanck())		

    def testVector01(self):
        dim = random.randint(1, 1000)		
        vec = H.Vector(dim)
        test = numpy.all(vec.data == numpy.zeros(dim))
        self.assertTrue(test)

    def testVector02(self):
        dim = random.randint(1,1000)
        x,y = self._set_xy(dim)
        vec = H.Vector(dim, x)
        test = numpy.all(vec.data == x )
        self.assertTrue(test)

    def testVector_norm(self):
        dim = random.randint(1,1000)    
        x,y = self._set_xy(dim)
        vec = H.Vector(dim,x)
        npnorm = numpy.sum(y*numpy.conj(y))
        norm = vec.norm()
        self.assertAlmostEqual(norm, npnorm.real)
    def testVector_inner(self):
        dim = random.randint(1,1000)    
        x, y= self._set_xy(dim)
        xx, yy= self._set_xy(dim)        
        vec = H.Vector(dim,x)
        minner = vec.inner(xx)
        ninner = numpy.sum(y*yy)
        self.assertAlmostEqual(minner.real, ninner.real)
        self.assertAlmostEqual(minner.imag, ninner.imag)        

    def _set_xy(self,dim):
        xre = numpy.random.random(dim) -0.5
        xim = numpy.random.random(dim) - 0.5
        x = [mpmath.mpc(xre[i],xim[i]) for i in range(dim)]
        y = xre + 1.j*xim
        return x, y


if __name__ == '__main__':
    unittest.main()
