import unittest
import random
import mpmath
import numpy
import time
import mpArray 


class TestQmap(unittest.TestCase):
    def setUp(self):
        mpmath.mp.dps = 40
        self.verbose = False
    
    def test_mpArray_init_(self):
        n = random.randint(1,5000)
        #x = [random.random() for i in range(n)]
        #self.assertRaises(TypeError, mpArray.mpArray, x)
        
        y = [mpmath.rand() for i in range(n)]
        arr = mpArray.mpArray(y)
        self.assertTrue(numpy.all(arr==y))
        
        z = [mpmath.mpc(mpmath.rand(),mpmath.rand()) for i in range(n)]
        
        z_real = [zz.real for zz in z]
        z_imag = [zz.imag for zz in z]
        arr = mpArray.mpArray(z)
        self.assertTrue(numpy.all(arr==z))
        self.assertTrue(numpy.all(arr.real() == z_real))
        self.assertTrue(numpy.all(arr.imag() == z_imag))
        
        start = time.clock()
        x=[mpmath.mpc('0','0')]*n
        t1 = time.clock() - start
        start = time.clock()
        x=[mpmath.mpc('0','0') for i in range(n)]
        t2 = time.clock() -start
        if self.verbose:
            print("[]*",n,t1,"vs.[for i in range(n)]",n,t2)
        
        
    def test_mpArray_add(self):
        n = random.randint(1,5000)
        x1 = [mpmath.mpc(mpmath.rand(), mpmath.rand()) for i in range(n)]
        x2 = [mpmath.mpc(mpmath.rand(), mpmath.rand()) for i in range(n)]        
        arr1 = mpArray.mpArray(x1)
        arr2 = mpArray.mpArray(x2)
        
        start = time.clock()
        x3 = [x1[i] + x2[i] for i in range(n)]
        t1 = time.clock() -start
        start = time.clock() 
        arr3 = arr1 + arr2
        t2 = time.clock() -start
        self.assertTrue(numpy.all(x3==arr3))
        if self.verbose:        
            print("-----------------------")
            print("number of ",n,"data add(list) time:", t1)
            print("number of ",n,"data add(mpArray) time:", t2)
    
    def test_mpArray_sub(self):
        n = random.randint(1,5000)
        x1 = [mpmath.mpc(mpmath.rand(), mpmath.rand()) for i in range(n)]
        x2 = [mpmath.mpc(mpmath.rand(), mpmath.rand()) for i in range(n)]        
        arr1 = mpArray.mpArray(x1)
        arr2 = mpArray.mpArray(x2)
        
        start = time.clock()
        x3 = [x1[i] - x2[i] for i in range(n)]
        t1 = time.clock() -start
        start = time.clock() 
        arr3 = arr1 - arr2
        t2 = time.clock() -start
        self.assertTrue(numpy.all(x3==arr3))
        if self.verbose:        
            print("-----------------------")        
            print("number of ",n,"data sub(list) time:", t1)
            print("number of ",n,"data sub(mpArray) time:", t2)
    
    def test_mpArray_mul(self):
        n = random.randint(1,5000)
        x1 = [mpmath.mpc(mpmath.rand(), mpmath.rand()) for i in range(n)]
        x2 = [mpmath.mpc(mpmath.rand(), mpmath.rand()) for i in range(n)]        
        arr1 = mpArray.mpArray(x1)
        arr2 = mpArray.mpArray(x2)
        
        start = time.clock()
        x3 = [x1[i] - x2[i] for i in range(n)]
        t1 = time.clock() - start
        start = time.clock() 
        arr3 = arr1 - arr2
        t2 = time.clock() -start
        self.assertTrue(numpy.all(x3==arr3))
        if self.verbose:        
            print("-----------------------")        
            print("number of ",n,"data mul(list) time:", t1)
            print("number of ",n,"data mul(mpArray) time:", t2)
    
    def test_mpArray_div(self):
        n = random.randint(1,5000)
        x1 = [mpmath.mpc(mpmath.rand(), mpmath.rand()) for i in range(n)]
        x2 = [mpmath.mpc(mpmath.rand(), mpmath.rand()) for i in range(n)]        
        arr1 = mpArray.mpArray(x1)
        arr2 = mpArray.mpArray(x2)
        
        start = time.clock()
        x3 = [x1[i] / x2[i] for i in range(n)]
        t1 = time.clock() - start
        start = time.clock() 
        arr3 = arr1 / arr2
        t2 = time.clock() -start
        self.assertTrue(numpy.all(x3==arr3))
        if self.verbose:        
            print("-----------------------")
            print("number of ",n,"data mul(list) time:", t1)
            print("number of ",n,"data mul(mpArray) time:", t2)
    
    def test_mpArray_index(self):
        n = random.randint(1,5000)
        x = mpmath.linspace(0,1,n)
        arr = mpArray.mpArray(x)
        start = time.clock()
        index1 = [xx < 0.5 for xx in x]
        t1 = time.clock() - start
        start = time.clock()
        index2 = arr < 0.5
        t2 = time.clock() - start
        self.assertTrue(numpy.all(index1==index2))
        if self.verbose:        
            print("-----------------------")        
            print("x < 0.5:","number of",n,"data serarch(list)", t1)
            print("x < 0.5:","number of",n,"data serarch(mpArray)", t1)
        
    def test_mpArray_toarray(self):
        n = random.randint(1,5000)
        input = numpy.random.random(n) + 1.j*numpy.random.random()
        y = [mpmath.mpc(x.real, x.imag) for x in input]
        arr = mpArray.mpArray(y)
        out = arr.toarray()
        self.assertTrue(numpy.all(out == input))
                
    def test_Matrix01(self):
        a = mpmath.mpc(mpmath.pi,"0")
        b = mpmath.mpc("0","1")
        m = mpArray.mpArray([[a,b],[b,a]])
        mat = mpArray.mpMatrix(m)
        evals, evecs = mat.eigen(False, 'qeispack', verbose=True)
        
        self.assertTrue(mpmath.fabs(mpmath.mpc(mpmath.pi,"1")  - evals[0]) < 1e-32)
        self.assertTrue(mpmath.fabs(mpmath.mpc(mpmath.pi,"-1")  - evals[1]) < 1e-32)
        self.assertTrue(numpy.all([mpmath.fabs(mpmath.mpf(1)/mpmath.sqrt("2") - mpmath.fabs(x)) < 1e-32 for x in evecs[0]]))
        self.assertTrue(numpy.all([mpmath.fabs(mpmath.mpf(1)/mpmath.sqrt("2") - mpmath.fabs(x)) < 1e-32 for x in evecs[1]]))
        
        
        self.assertTrue(mpmath.fabs(mpmath.fsum(evecs[0]*evecs[0].conj())) - mpmath.mpf(1) < 1e-32)
        self.assertTrue(numpy.all([vec.abs2() -mpmath.mpf(1) < 1e-32 for vec in evecs]))
        self.assertTrue(mpmath.fabs(evecs[0].inner(evecs[1])) < 1e-32)

 