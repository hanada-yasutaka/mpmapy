import unittest
import random
import mpmath
import numpy
import time
import mpmapy
from mpmapy import QmapSystem
from mpmapy.QmapSystem import ScaleInfo
from mpmapy import maps

#import QmapSystem as QS
#import mpMaps as Maps

class QmapTest(unittest.TestCase):
    def setUp(self):
        mpmath.mp.dps = 40
        k=mpmath.rand()*mpmath.mpf("10")
        self.map = maps.StandardMap(k=k)
        
    def testScaleInfo01(self):
        dim = 4
        domain = [[0,1],[0,1]]
        scl = ScaleInfo(dim, domain)
        a = (mpmath.mpf(0), mpmath.mpf(1))
        self.assertTrue(scl.domain == [a,a]) 
        self.assertTrue( numpy.all(scl.x[0] == mpmath.linspace(0,1, dim, endpoint=False)) )
        self.assertTrue( numpy.all(scl.x[1] == mpmath.linspace(0,1, dim, endpoint=False)) )
        self.assertTrue(scl.h == mpmath.mpf("0.25"))

    def testScaleInofo02(self):
        dim = 4
        domain = [["1","-1"],["0","1"]]
        self.assertRaises(AssertionError, ScaleInfo, dim, domain)
        domain = [["1","1"],["0","1"]]
        self.assertRaises(AssertionError, ScaleInfo, dim, domain)        
        domain = [[1,-1],["0","1"]]
        self.assertRaises(AssertionError, ScaleInfo, dim, domain)
        domain = [[1,1],["0","1"]]
        self.assertRaises(AssertionError, ScaleInfo, dim, domain)        
        domain = [[mpmath.mpf("1"),mpmath.mpf("1")],["0","1"]]
        self.assertRaises(AssertionError, ScaleInfo, dim, domain)                
        domain = [[mpmath.mpf("1"),mpmath.mpf("-1")],["0","1"]]
        self.assertRaises(AssertionError, ScaleInfo, dim, domain)                        
        
    def testQmapSystem(self):
        dim = random.randint(1,20)
        print dim
        qmin = -mpmath.rand()
        qmax = mpmath.rand()
        pmin = -mpmath.rand()
        pmax = mpmath.rand()
        domain = [[qmin,qmax],[pmin,pmax]]
        qmapsys = QmapSystem(map=self.map, type='U', dim=dim, domain=domain)
        

        #self.assertTrue( numpy.all(qmapsys.qmap.scaleinfo.x[0] == numpy.linspace(qmin,qmax, dim, endpoint=False)) )
        #self.assertTrue( numpy.all(qmapsys.qmap.scaleinfo.x[1] == numpy.linspace(pmin,pmax, dim, endpoint=False)) )
    
    def __test_getEigen(self):
        dim = random.randint(1,20)

        qmin = -mpmath.rand()
        qmax = mpmath.rand()
        pmin = -mpmath.rand()
        pmax = mpmath.rand()
        domain = [[qmin,qmax],[pmin,pmax]]
        qmapsys = QmapSystem(map=self.map, type='U', dim=dim, domain=domain)
        evals, evecs = qmapsys.getEigen()
        for evec in evecs:
            self.assertTrue(mpmath.fabs(evec.norm() -mpmath.mpf(1) ) < 1e-32)
            
        for i in range(len(evecs)):
            vals = [ evecs[i].inner(evecs[j]) for j in range(len(evecs)) ]
            self.assertTrue(mpmath.fabs(mpmath.fabs(vals[i]) - 1.0) < 1e-32) 
            vals.pop(i)
            index = [x < 1e-32 for x in vals]
            self.assertFalse(numpy.all(index))
        for evec in evecs:
            evec.hsmrep(10,10)
