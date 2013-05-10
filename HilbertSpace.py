# -*- coding: utf-8 -*-
import mpmath
import mpfft
import time

def sort_nicely( l ): 
	""" 
    Sort the given list in the way that humans expect. For example,
    
    >>> l = ["1", "4", "2" , "0"]
    >>> sort_nicely(l)
    >>> l
    ['0', '1', '2', '4']
    
    Note that negative values (its argument takes only string!) is shifted to the backward  

    >>> l = ["-1","3","-2", "4"]
    >>> mapy.QmapSystem.sort_nicely(l)
    >>> l
    ['3', '4', '-1', '-2']
     
    """ 
	import re     
	convert = lambda text: int(text) if text.isdigit() else text 
	alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
	l.sort( key=alphanum_key )    



class HilbertSpace(object):
	def __init__(self, dim):
		if not isinstance(dim, int):
			raise ValueError('Hilbert dim has to be an integer value')
		self.dim = dim

class _Vector(HilbertSpace):
	def __init__(self, dim, data):
		if isinstance(data, _Vector):
			self._data = data._data
		elif isinstance(data, list):
			self._data = [mpmath.mpf(x) for x in data]
		elif data == None and isinstance(dim, int):
			self._data = [mpmath.mpf(0)] * dim
		else:
			raise ValueError("")						
		HilbertSpace.__init__(self, len(self._data))
		self.data = getattr(self, '_data')
	def __repr__(self):
		return repr(self._data)
		
	def __str__(self):
		return str(self._data)
		
	def __call__(self):
		return self._data
	
	def __add__(self, y):
		return self.add(y)
		
	def __sub__(self, y):
		return self.sub(y)

	def __mul__(self, y):
		return self.mul(y)
		
	def __div__(self, y):
		return self.div(y)
			
	def __len__(self):
		return len(self._data)
	def add(self, y):
		if not isinstance(y, Vector):
			y = _Vector(self.dim, y)
		res = [mpmath.fadd(self._data[i], y._data[i]) for i in range(self.dim)]
		return _Vector(self.dim, res)

	def sub(self, y):
		if not isinstance(y, Vector):
			y = _Vector(self.dim, y)
		res = [mpmath.fsub(self._data[i], y._data[i]) for i in range(self.dim)]
		return _Vector(self.dim, res)
	def mul(self, y):
		if not isinstance(y, Vector):
			y = _Vector(self.dim, y)
		res = [mpmath.fmul(self._data[i], y._data[i]) for i in range(self.dim)]
		return _Vector(self.dim, res)
	def div(self, y):
		# todo zero divided
		if not isinstance(y, Vector):
			y = _Vector(self.dim, y)
		res = [mpmath.fdiv(self._data[i], y._data[i]) for i in range(self.dim)]
		return _Vector(self.dim, res)
		
class Vector(_Vector):
	def __init__(self, input, imag_data=None):
		#_Vector.__init__(self, dim, None)
		if isinstance(input, Vector):
			self._data = input._data
		elif isinstance(input, int):
			self._data = [mpmath.mpc("0","0")] * input
		elif isinstance(input, list):
			if (imag_data == None):
				self._data = [mpmath.mpc(x.real, x.imag) for x in input]
			elif len(input) != len(imag_data):
				raise ValueError("input data must be same length.")				
			else: #(imag_data != None):
				self._data =[mpmath.mpc(input[i], imag_data[i]) for i in range(len(input))]
		else:
			raise ValueError("")
		self.dim = len(self._data)

		self.data = getattr(self, '_data')
		
		
	def add(self, y):
		if not isinstance(y, Vector):
			y = Vector(y)
		res = [mpmath.fadd(self._data[i], y._data[i]) for i in range(self.dim)]
		return Vector(res)

	def sub(self, y):
		if not isinstance(y, Vector):
			y = Vector(y)
		res = [mpmath.fsub(self._data[i], y._data[i]) for i in range(self.dim)]
		return Vector(res)
	def mul(self, y,test=False):

		if not isinstance(y, Vector):
			y = Vector(y)
		res = [mpmath.fmul(self._data[i] , y._data[i]) for i in range(self.dim)]
		return Vector(res)
	
	def div(self, y):
		# todo zero divided
		if not isinstance(y, Vector):
			y = Vector(y)
		res = [mpmath.fdiv(self._data[i], y._data[i]) for i in range(self.dim)]
		return Vector(res)
		
	def one(self,i):
		self._data[i] = mpmath.mpc("1","0")
	
	def set(self,i,x):
		if not isinstance(x,mpmath.mpc) and not isinstance(x, mpmath.mpf):
			raise ValueError("excepted mpmath.mpc or mpmath.mpf value")
		self._data[i] = x
	
	def setall(self,x):
		self._data = [self.set(i,x[i]) for i in range(len(x)) ] 
	
	def val(self):
		return self._data		
		
	def norm(self):
		norm = mpmath.fsum([mpmath.fmul(x, mpmath.conj(x)) for x in self._data])
		return norm.real	

	def conj(self):
		return [mpmath.conj(x) for x in self._data]
		
	def dot(self, vec):
		if isinstance(vec, Vector):
			return mpmath.fsum(mpmath.fmul(self._data[i], vec.data[i]) for i in range(self.dim))
		else:
			return mpmath.fsum(mpmath.fmul(self._data[i], vec[i]) for i in range(self.dim))

	def inner(self,vec):
		if isinstance(vec, Vector):
			return mpmath.fsum([self._data[i]*mpmath.conj(vec.data[i]) for i in range(self.dim)])
		else:
			return mpmath.fsum([self._data[i]*mpmath.conj(vec[i]) for i in range(self.dim)])

	def abs2(self):
		return [mpmath.fabs(mpmath.fmul(self._data[i],mpmath.conj(self._data[i]))) for i in range(self.dim)]
						
	def real(self):
		res = [ x.real for x in self._data ] 
		return _Vector(self.dim, res)

	def imag(self):
		res = [ x.imag for x in self._data ]
		return _Vector(self.dim, res)


class Matrix(HilbertSpace):
	def __init__(self, dim, data=None):
		HilbertSpace.__init__(self,dim)  
	
		if data != None:
			if not isinstance(data, mpmath.matrix):
				raise TypeError('matrix type must be mpmath.matrix')
			if not (data.cols, data.rows) == (dim, dim):
				raise ValueError('matrix dimension must be (%d,%d)' % (dim, dim))

		self.data = data if data != None else mpmath.zeros(self.dim, self.dim)
		self.evals = None
		self.evecs = None
	
	def eigen(self,left=False, solver='qeispack', verbose=False):
		solver_list = {'qeispack':self.qeispack, 'lapack':self.lapack}
		if solver not in solver_list:
			raise TypeError("excepted solver:", solver_list.keys())
		self.solver = solver_list[solver]
		self.solver(verbose)

	def qeispack(self,verbose):
		from utility import ctypes_wrapper as cw
	
		start=time.time()
		rfname = 'matrix_real.dat'
		ifname = 'matrix_imag.dat'
		self._save_matrix(rfname, ifname)
		end = time.time()
		t = time.time() - start
		if verbose:
			print("save matrix size", self.dim, "in", t,'sec.')			

		start= time.time()
		qeispack = cw.qEispack()
		qeispack.file2call_eig(self.dim, rfname, ifname)
		end = time.time()
		t = time.time() -start
		
		self._load_file(verbose)
		if verbose:
			print("Computation eigen-values and -vectors of size", self.dim, "in", t,'sec.')			

	def lapack(self,verbose):
		pass

	def _save_matrix(self, rfname, ifname):
		file1 = open(rfname, "w")
		file2 = open(ifname, "w")
		for i in range(self.data.cols):
			[file1.write("%s " % self.data[i,j].real) for j in range(self.data.rows)]
			[file2.write("%s " % self.data[i,j].imag) for j in range(self.data.rows)]		
			file1.write("\n")
			file2.write("\n")		
		file1.close()
		file2.close()
		
	def _load_file(self,verbose=False):
		import glob, os
		from utility import sort
		list = glob.glob('eigen_vec*.dat')
		sort.natural(list)
		self.evals = [None for i in range(self.dim)]
		self.evecs = [None for i in range(self.dim)]
		start = time.time()
		for i,fname in enumerate(list):
			f = open(fname, "r")
			data = [mpmath.mpc(x.split(",")[0],x.split(",")[1]) for x in f]
			self.evecs[i] = Vector(data) 
			f.close()
			os.remove(fname)
		f = open("eigen_vals.dat")
		data = [ mpmath.mpc(x.split(",")[0], x.split(",")[1]) for x in f ]
		self.evals = Vector(data)
		f.close()
		os.remove("eigen_vals.dat")
		t = time.time() - start
		if verbose:
			print("load %d eigen-value and -vector:" ,t, "sec.")

	def getEigen(self):
		return self.evals, self.evecs

