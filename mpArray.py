import numpy
import mpmath
import time
class mpArray(numpy.ndarray):
    def __new__(cls, input, dtype='object'):
        if isinstance(input, int):
            data = [mpmath.mpc('0','0')]*input
            obj = numpy.asarray(data, dtype='object').view(cls)
        else:
            obj = numpy.asarray(input, dtype=dtype).view(cls)
        return obj
    def __array_finalize__(self, obj):
        if obj is None: return
    
    def real(self):
        #return mpArray([ x.real for x in self])
        return numpy.asarray([x.real for x in self], dtype='object').view(type(self))

    def imag(self):
        return numpy.asarray([x.imag for x in self], dtype='object').view(type(self))
        #return mpArray([ x.imag for x in self])
    
    def toarray(self):
        return numpy.array(self.tolist(),dtype=numpy.complex128)
    
    def norm(self):
        return mpmath.fabs(self.inner())
    
    def inner(self, out_arr=None):
        if out_arr !=None and not isinstance(out_arr, mpArray):
            raise ValueError("input data must be mpArray") 
        if out_arr==None:
            return  self.dot(numpy.conj(self))
        else:
            return self.dot(numpy.conj(out_arr))
    def abs2(self):
        return mpArray(numpy.abs(self*numpy.conj(self)))



class mpMatrix(numpy.ndarray):
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
        
        return self._load_file(verbose)
        if verbose:
            print("Computation eigen-values and -vectors of size", len(self), "in", t,'sec.')            

    def lapack(self,verbose):
        matrix = self.toarray()

        evals, evecs = numpy.linalg.eig(matrix)
        return evals, evecs.transpose()
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
    
