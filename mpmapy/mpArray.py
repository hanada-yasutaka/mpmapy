# -*- coding: utf-8 -*-
import numpy
import mpmath
import time
class mpArray(numpy.ndarray):
    """
    
    Description
    ==========
    
    mpmath.mpcをnumpy.ndarrayでラップする為のクラスです
    
    Create new numpy.ndarray
    
    
    .. math:: (a + b)^2 = a^2 + 2ab + b^2
    
    
    Parameters
    ----------
    
    input : int or list or array
        if input is integer, new array filled zeros, 
        else new array filled input values.
    dtype : date-type, (optional)
    
    
    Returns
    -----------
    
    our : ndarray
        An array object satisfying the specified requirements
    
    Attributes
    ----------
    real : mpArray
        Real part of the array
    imag : mpArray
        Imaginary part of the array
    
    See Also
    ----------
    
    numpy.ndarray, numpy.ndarray subclassing <http://docs.scipy.org/doc/numpy/user/basics.subclassing.html>
    
    .. note::
    
        echo "test"
    
    Examples
    ---------
    
    >>> import mpmapy
    >>> print (mpmapy.mpArray(4))
    [(0.0 + 0.0j) (0.0 + 0.0j) (0.0 + 0.0j) (0.0 + 0.0j)]
    >>> print (mpmapy.mpArray([1,2,3]))
    [1 2 3]
    >>> y = mpmapy.mpArray.linspace(0,1,4) + mpmapy.mpArray.linspace(1,2,4)*1.j
    >>> y = mpmapy.mpArray.linspace(0,1,4) + mpmapy.mpArray.linspace(1,2,4)*1.j
    >>> print y.real
    [0.0 0.25 0.5 0.75]
    >>> print(y.real)
    [0.0 0.25 0.5 0.75]
    >>> print(y.imag)
    [1.0 1.25 1.5 1.75]
    >>> y.toarray(), y.toarray().dtype
    (array([ 0.00+1.j  ,  0.25+1.25j,  0.50+1.5j ,  0.75+1.75j]), dtype('complex128'))
    
    >>> # todo
    Traceback (most recent call last):
      File "/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/doctest.py", line 1289, in __run
        compileflags, 1) in test.globs
      File "<doctest __main__.mpArray[9]>", line 1, in <module>
        y = mpmapy.mpArray.linspace(0,1,4) + 1.j*mpmapy.mpArray.linspace(1,2,4)
    TypeError: unsupported operand type(s) for *: 'complex' and 'mpArray'
    
    """

    def __new__(cls, input=[], dtype='object'):
        if isinstance(input, int):
            data = [mpmath.mpc('0','0')]*input
            obj = numpy.asarray(data, dtype='object').view(cls)
        else:
            obj = numpy.asarray(input, dtype=dtype).view(cls)
        return obj

    def __array_finalize__(self, obj):
        if obj is None: return
    
    @property
    def real(self):
        return numpy.asarray([x.real for x in self], dtype='object').view(type(self))
    
    @property
    def imag(self):
        return numpy.asarray([x.imag for x in self], dtype='object').view(type(self))
    
    def toarray(self):
    	# to numpy array, dtype is numpy.complex128
        return numpy.array(self.tolist(),dtype=numpy.complex128)
    
    def insert(self, i, x):
        """
        >>> import mpmapy
        >>> import mpmath
        >>> x = mpmapy.mpArray(2)
        >>> x.insert(0,mpmath.mpc("1","1"))
        >>> print(x)
        [(1.0 + 1.0j) (0.0 + 0.0j)]
        """
        self[i] = x

    def norm(self):
        """
        >>> import mpmapy
        >>> import mpmath
        >>> x = mpmapy.mpArray([1, -1])
        >>> print x.norm()
        1.4142135623731
        >>> print x.norm() == mpmath.sqrt("2")
        True
        """
        return mpmath.sqrt(mpmath.fabs(self.inner()))
    def norm2(self):
        return mpmath.fabs(self.inner())
    
    def inner(self, out_arr=None):
        """
        
        """
        if out_arr !=None and not isinstance(out_arr, mpArray):
            raise ValueError("input data must be mpArray")
        if out_arr==None:
            return  self.dot(numpy.conj(self))
        else:
            return self.dot(numpy.conj(out_arr))
    
    def fftshift(self):
        return mpArray(numpy.fft.fftshift(self))
    
    def abs2(self):
        return mpArray(numpy.abs(self*numpy.conj(self)))
    
    def normalize(self):
        """
        >>> import mpmapy
        >>> import mpmath
        >>> x = mpmapy.mpArray([1, -1])
        >>> print(x.normalize())
        [0.707106781186547 -0.707106781186547]
        """
        return self/self.norm()

    def conj(self):
        return numpy.conj(self)
    
    def parity(self):
        vec = self.tolist()
        vec.append(self[0])
        vec1 = mpArray(vec)
        vec1 = vec1.normalize()#/mpmath.sqrt(vec1.norm())
        vec2 = mpArray([x for x in vec1[::-1]])
        return vec1.inner(vec2).real

    @classmethod
    def ones(cls, dim):
        return mpArray([mpmath.mpc("1", "0")]*dim)
    @classmethod
    def zeros(cls, dim):
        return mpArray(dim)
    @classmethod
    def linspace(cls, min, max, num, endpoint=False):
        return mpArray(mpmath.linspace(min, max, num, endpoint=endpoint))

        
class mpMatrix(numpy.ndarray):
    """
    Examples
    ----------
    >>> import mpmapy
    >>> import mpmath
    >>> mpmath.mp.dps=36
    >>> print(mpmapy.mpMatrix(2))
    [[(0.0 + 0.0j) (0.0 + 0.0j)]
     [(0.0 + 0.0j) (0.0 + 0.0j)]]
    >>> mat= mpmapy.mpMatrix([[mpmath.mpc('0','1'), mpmath.pi],[mpmath.pi,mpmath.mpc('0','1')]])
    >>> print(mat)
    [[(0.0 + 1.0j) 3.14159265358979323846264338327950288]
     [3.14159265358979323846264338327950288 (0.0 + 1.0j)]]
    >>> eval, evec = mat.eigen()
    >>> print(eval)
    [(3.1415926535897932384626433832795028 + 1.0j)
     (-3.1415926535897932384626433832795028 + 1.0j)]
    >>> print(evec)
    [mpArray([(0.707106781186547524400844362104848992 + 0.0j),
           (0.707106781186547524400844362104848992 + 0.0j)], dtype=object), mpArray([(-0.707106781186547524400844362104848992 + 0.0j),
           (0.707106781186547524400844362104848992 + 0.0j)], dtype=object)]
    >>> eval[0] - mpmath.mpc(mpmath.pi,"1")
    mpc(real='-8.4259435066941568571199034988090570186e-35', imag='0.0')
    """
    def __new__(cls, input, dtype='object'):
        if isinstance(input, int):
            data  = [[mpmath.mpc("0","0")]*input]*input
            obj = numpy.asarray(data, dtype=object).view(cls)
        else:
            obj = numpy.asarray(input, dtype=dtype).view(cls)            
            if len(obj.shape) != 2 and obj.shape[0] != obj.shape[1]:
                raise TypeError("excepted: n by n matrix, but input data shape",obj.shape )
        return obj
    
    def eigen(self,left=False, solver='qeispack', verbose=False):
        solver_list = {'qeispack':self.qeispack, 'lapack':self.lapack}
        if solver not in solver_list:
            raise TypeError("excepted solver:", solver_list.keys())
        self.solver = solver_list[solver]
        return self.solver(verbose)
        

    def qeispack(self,verbose):
        from utility import ctypes_wrapper as cw

        start=time.time()
        rfname = 'matrix_real.dat'
        ifname = 'matrix_imag.dat'
        self._save_matrix(rfname, ifname)
        end = time.time()
        t = time.time() - start
        if verbose:
            print("save matrix size", self.shape, "in", t,'sec.')            

        start= time.time()
        qeispack = cw.qEispack()
        qeispack.file2call_eig(len(self), rfname, ifname)
        end = time.time()
        t = time.time() -start
        if verbose:
            print("(QEISPACK) get eigen-values and -vectors of size", len(self), "in", t,'sec.')
        return self._load_file(verbose)


    def lapack(self,verbose):
        start = time.time()        
        matrix = self.toarray()

        evals, evecs = numpy.linalg.eig(matrix)
        evals = mpArray(evals)
        evecs = [mpArray(evec) for evec in evecs.transpose()]
        if verbose:
            t = time.time() -start
            print("(LAPACK) get eigen value and eigen vector",len(self), "in",t,'sec.')
        return evals, evecs 
    def toarray(self):
        return numpy.array(self.tolist(), dtype=numpy.complex128)
    
    def _save_matrix(self, rfname, ifname):
        file1 = open(rfname, "w")
        file2 = open(ifname, "w")
        for i in range(self.shape[0]):
            [file1.write("%s " % self[i,j].real) for j in range(self.shape[1])]
            [file2.write("%s " % self[i,j].imag) for j in range(self.shape[1])]        
            file1.write("\n")
            file2.write("\n")        
        file1.close()
        file2.close()
        
    def _load_file(self,verbose=False):
        import glob, os
        from utility import sort
        list = glob.glob('eigen_vec*.dat')
        sort.natural(list)
        evals, evecs  = [None] * len(self), [None]*len(self)

        start = time.time()

        for i,fname in enumerate(list):
            f = open(fname, "r")
            data = [mpmath.mpc(x.split(",")[0],x.split(",")[1]) for x in f]
            evecs[i] = mpArray(data) #Vector(data) 
            f.close()
            os.remove(fname)
        f = open("eigen_vals.dat")
        data = [ mpmath.mpc(x.split(",")[0], x.split(",")[1]) for x in f ]
        evals = mpArray(data)# Vector(data)
        f.close()
        os.remove("eigen_vals.dat")
        t = time.time() - start
        if verbose:
            print("load ",self.shape," eigen-value and -vector:" ,t, "sec.")
        return evals, evecs
    
 

def _test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    _test()
