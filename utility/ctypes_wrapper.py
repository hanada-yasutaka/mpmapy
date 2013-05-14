import ctypes
import os
import numpy

class qEispack(object):
	def __init__(self):
		path = os.environ['PYTHONPATH'].split(":")
		index = [p.endswith("mpmapy") for p in path].index(True)
		file_path = path[index] + '/shared/wrapper_qeispack.so'
		self.module = ctypes.cdll.LoadLibrary(file_path)		
	
	def file2call_eig(self, dim, rfname, ifname	):
		self.module.file2call_eig_.argtypes = [ ctypes.POINTER(ctypes.c_int),\
												ctypes.POINTER(ctypes.c_int),\
												ctypes.c_char_p, \
												ctypes.c_char_p ]

		self.module.file2call_eig_.restype = ctypes.c_void_p
		length1=len(rfname)
		length2=len(rfname)

		if length1 != length2:
			raise ValueError
		# todo pytthon3 TypeError: bytes or integer address expected instead of str instance
		rfname = ctypes.c_char_p(rfname)
		ifname = ctypes.c_char_p(ifname)

		dim = ctypes.c_int(dim)
		length=ctypes.c_int(length1)
		self.module.file2call_eig_(ctypes.byref(dim), ctypes.byref(length),rfname, ifname)




class call_hsm_rep(object): 
	def __init__(self):
		path = os.environ['PYTHONPATH'].split(":")
		index = [p.endswith('mpmapy') for p in path].index(True)
		file_path = path[index] + '/shared/libhsm.so'
		self.c_lib = ctypes.cdll.LoadLibrary(file_path)
		
	def husimi_rep(self, vec, dim, domain, hsm_range, hsm_grid):
		terget = numpy.array(vec,dtype=numpy.complex128)
		qmin, qmax = domain[0][0], domain[0][1]
		pmin, pmax = domain[1][0], domain[1][1]
		#dim= dim
		h = (qmax - qmin)*(pmax - pmin)/dim
		vqmin, vqmax = hsm_range[0][0], hsm_range[0][1]
		vpmin, vpmax = hsm_range[1][0], hsm_range[1][1]
		row, col = hsm_grid[0], hsm_grid[1]
		hsm_data = numpy.zeros([row, col], dtype=numpy.float64)
		
		rvec = numpy.copy(vec.real)
		ivec = numpy.copy(vec.imag)
		
		f4ptr = ctypes.POINTER(ctypes.c_double)
		data = (f4ptr*row)(*[ROW.ctypes.data_as(f4ptr) for ROW in hsm_data])
		
		hsm = self.c_lib.wrapper_husimi_rep
		hsm.restype = None
		hsm.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),\
                        ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.c_int, 
                        ctypes.c_int, ctypes.c_int,\
                        ctypes.c_double, ctypes.c_double,\
                        ctypes.c_double, ctypes.c_double,\
                        ctypes.c_double,\
                        ctypes.c_double, ctypes.c_double,\
                        ctypes.c_double, ctypes.c_double]
		hsm(rvec.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),  \
             ivec.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),  \
             data,\
             len(terget),\
             row, col,\
             qmax, qmin,\
             pmax, pmin,\
             h,\
             vqmax, vqmin,\
             vpmax, vpmin
            )
		hsm_data=numpy.array([data[i][j] for i in range(row) for j in range(col)])
		hsm_imag = hsm_data.reshape(row,col)
#		x = numpy.linspace(vqmin, vqmax, row)
#		y = numpy.linspace(vpmin, vpmax, col)
#		X,Y = numpy.meshgrid(x,y)

		return hsm_imag
