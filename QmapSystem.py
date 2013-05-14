from mpArray import *
import mpmath
import mpfft
twopi = mpmath.mpf("2")*mpmath.pi

class HilbertSpace(object):
    def __init__(self, dim):
        if not isinstance(dim, int):
            raise ValueError('excepted integer')
        self.dim = dim

class ScaleInfo(HilbertSpace):
    def __init__(self, dim):
        HilbertSpace.__init__(self, dim)
        self.domain = [[None,None] for i in range(2)]
        
    def setQdomain(self, qmin, qmax):
        if not isinstance(qmin,int) and not isinstance(qmin, mpmath.mpf):
            raise ValueError('input values have to be an integer or an mpmath.mpf value')
        if not isinstance(qmax,int) and not isinstance(qmax, mpmath.mpf):
            raise ValueError('input values have to be an integer or an mpmath.mpf value')
        if qmin >= qmax:
            raise ValueError('Error: qmin >= qmax ')
        self.domain[0][0] = mpmath.mpf(qmin)
        self.domain[0][1] = mpmath.mpf(qmax)
    
    def setPdomain(self, pmin, pmax):
        if not isinstance(pmin,int) and not isinstance(pmin, mpmath.mpf):
            raise ValueError('input values have to be an integer or an mpmath.mpf value')
        if not isinstance(pmax,int) and not isinstance(pmax, mpmath.mpf):
            raise ValueError('input values have to be an integer or an mpmath.mpf value')
        if pmin > pmax:
            raise ValueError('Error: pmin >= pmax')
        self.domain[1][0] = mpmath.mpf(pmin)
        self.domain[1][1] = mpmath.mpf(pmax)

    def setDomain(self, r):
        self.setQdomain(r[0][0], r[0][1])
        self.setPdomain(r[1][0], r[1][1])
        self.x = [mpmath.linspace(self.domain[i][0], self.domain[i][1], self.dim, endpoint=False) for i in range(2)]
        
    def getDomain(self):
        return self.domain

    def getX(self):
        return self.x

    def getVolume(self):
        return (self.domain[0][1] - self.domain[0][0])* (self.domain[1][1] - self.domain[1][0])
    
    def getPlanck(self):
        return self.getVolume()/self.dim

class State(mpArray):
    def __new__(cls, scaleinfo, data=None):
        if data != None and isinstance(data, mpArray):
            obj = numpy.asarray(data, dtype='object').view(cls)
        elif isinstance(scaleinfo, ScaleInfo):
            data = [mpmath.mpc('0','0')]*scaleinfo.dim
            obj = numpy.asarray(data, dtype='object').view(cls)
        else:
            raise TypeError("usage: scaleinfo, (optional) mpArray")
        if len(obj) != scaleinfo.dim:
            raise ValueError("data length does not match")            
        obj.scaleinfo = scaleinfo
        obj.dim = scaleinfo.dim
        return obj
    def __array_finalize__(self, obj):
        if obj is None: return
        self.scaleinfo = getattr(obj, 'scaleinfo',None)
        self.dim = getattr(obj, 'dim', None)
        
    def qrep(self):
        return self
    
    def prep(self):
        pass
    
    def parity(self):
        vec = [x for x in self]
        vec.append(self[0])
        arr = mpArray(vec)
        return arr.inner(arr[::-1]).real 
    
    def set_qconst(self, q_c):
        self = State(self.scaleinfo)
        q = self.scaleinfo.x[0]
        dx = mpmath.fabs(q[1] - q[0])
        qmin = self.scaleinfo.domain[0][0]
        qmax = self.scaleinfo.domain[0][1]
        index = int(round((q_c - qmin)/dx))
        if index > self.dim - 1 or index < 0:
            raise ValueError('q_c set in (%s,%s) with %s' % (mpmath.nstr(qmin), mpmath.nstr(qmax), mpmath.nstr(dx)) )
        self[index] = mpmath.mpc("1", "0")
    
    def hsmrep(self, col, row, hsm_region=None):
        if hsm_region==None:
            self._hsm_region = self.scaleinfo.domain
        else:
            self._hsm_region = hsm_region 
        from utility import ctypes_wrapper
        cw = ctypes_wrapper.call_hsm_rep()
        hsm_imag = cw.husimi_rep(self, self.dim, self.scaleinfo.domain, self._hsm_region, [row,col])
        x = numpy.linspace(self._hsm_region[0][0], self._hsm_region[0][1], row)
        y = numpy.linspace(self._hsm_region[1][0], self._hsm_region[1][1], col)
        X,Y = numpy.meshgrid(x,y)       
        return X,Y,hsm_imag
    
    def _mult_save_hsmrep(self, **kwargs):
        col = kwargs['col']
        row = kwargs['row']
        title = kwargs['title']
        if 'region' in kwargs:
            hsm_region = kwargs['region']
        elif 'hsm_region' in kwargs:
            hsm_region = kwargs['hsm_region']
        else:
            hsm_region = None        
        self.save_hsmrep(col=col, row=row,hsm_region=hsm_region,title=title)
    
    def save_hsmrep(self,col, row, hsm_region=None, title=None):
        if title==None:
            raise ValueError("title is None")
        domain  = self.scaleinfo.domain
        import datetime   
        X, Y, hsm_imag = self.hsmrep(col,row, hsm_region)
        with open(title, "w") as of:
            of.write("# DATE %s\n" % datetime.datetime.now())
            of.write("# DIM %d\n# QMIN %f\n# QMAX %f\n" % (self.dim, domain[0][0], domain[0][1]))
            of.write("# PMIN %f\n# PMAX %f\n" % (domain[1][0],domain[1][1]))
            of.write("# VQMIN %f\n# VQMAX %f\n# VPMIN %f\n# VPMAX %f\n" % (self._hsm_region[0][0], self._hsm_region[0][1], self._hsm_region[1][0], self._hsm_region[1][1]))
            of.write("# ROW %d\n# COL %d\n" % (row, col))
            for slice_data in hsm_imag.transpose():
                numpy.savetxt(of, slice_data)
            of.write("\n") 


class Qmap(HilbertSpace):
    def __init__(self,dim):
        HilbertSpace.__init__(self, dim)
        self.scaleinfo = ScaleInfo(dim)
        self.stateIn = mpArray(dim)
        self.stateOut = mpArray(dim)
        
    def setIn(self, state):
        if len(state) != self.dim:
            if not isinstance(state, mpArray) or not isinstance(state, State):
                raise ValueError("expected State or mpArray or State (dim=%d) instance" % self.dim)
        self.stateIn = state

    def getIn(self):
        return State(self.scaleinfo, self.stateIn)
    
    def getOut(self):
        return State(self.scaleinfo, self.stateOut)

    def swap(self):
        tmp =self.stateIn.copy() #[self.stateIn.data[i] for i in range(self.dim)]
        self.stateIn  = self.stateOut.copy()
        self.stateOut = tmp

    def pull(self):
        self.stateIn = mpArray(self.stateOut)

    def operate(self):
        pass
    
    def map(self, state):
        self.setIn(state)
        self.operate()
        self.getOut()
    
    def evol(self):
        self.operate()
        self.pull()
        
class Unitary(Qmap):
    def __init__(self, map, dim, tau=1):
        Qmap.__init__(self, dim)
        self.map = map
        if not isinstance(tau, mpmath.mpf) and  not isinstance(tau, int):
            raise TypeError("tau types must be mpmath.mpf or mpmath.mpc")
    
        self.tau = tau
        self.matrix = None
        #self.evals = None
        #self.evecs = None
        
    def setDomain(self, r):
        self.scaleinfo.setDomain(r)
        self.h = self.scaleinfo.getPlanck()
        
    def op0(self, x, isShift=False):
        if isShift:
            sub_x = numpy.fft.fftshift(x)
            func = self.map.ifunc0(sub_x)
            self.operator[0] = mpArray([ mpmath.exp(mpmath.mpc("0", -twopi*func[i]*self.tau/self.h)) for i in range(self.dim) ])
        else:
            func = self.map.ifunc0(x)
            self.operator[0] = mpArray([ mpmath.exp(mpmath.mpc("0", -twopi*func[i]*self.tau/self.h)) for i in range(self.dim) ])

    def op1(self, x, isShift=True):
        if isShift:
            sub_x = numpy.fft.fftshift(x)
            func = self.map.ifunc1(sub_x)
            self.operator[1] = mpArray([ mpmath.exp(mpmath.mpc("0", -twopi*func[i]*self.tau/self.h)) for i in range(self.dim) ])
        else:
            func = self.map.ifunc1(x)
            self.operator[1] = mpArray([ mpmath.exp(mpmath.mpc("0", -twopi*func[i]*self.tau/self.h)) for i in range(self.dim) ])
        
    def operate(self, invec=None, queue=None, verbose=False):
        if invec != None:
            if isinstance(invec, mpArray):
                self.setIn(invec)
            elif isinstance(invec, list) and self.dim==len(invec):
                self.setIn(mpArray(invec))

        qvec = self.operator[0]*self.stateIn
        pvec = mpArray(mpfft.fft(qvec, inverse=False, verbose=verbose))
        pvec = self.operator[1]*pvec 
        self.stateOut = mpArray(mpfft.fft(pvec, inverse=True, verbose=verbose))
        if queue !=None: 
            queue.put(self.stateOut)
    
    def setOperator(self, fft=True):
        self.operator = [mpArray(self.dim) for i in range(2)] 
        
        self.op0(self.scaleinfo.x[0], False)
        self.op1(self.scaleinfo.x[1], self._isShift())        

    def setMatrix(self, fft=True, core=2, verbose=False):
        from utility import parallel
        self.setOperator(fft)
        matrix = mpMatrix(self.dim)
        

        orth_basis = [mpArray(self.dim) for i in range(self.dim)]
        
        for i in range(self.dim):
            orth_basis[i][i] = mpmath.mpc("1","0")

        start = time.time()
        
        if core > 1:
            #data = [orth_basis[i] for i in range(self.dim)]
            multi = parallel.MultiProcessing(self.operate, orth_basis, core=core)
            multi.run()
            res = multi.getResults()
            for i in range(self.dim):
                matrix[i,:] = res[i]
        elif core == 1:
            for i in range(self.dim):
                self.setIn(orth_basis[i])
                self.operate()
                
                matrix[i,:] = self.stateOut[:] #mpmath.matrix(self.stateOut.data).T
        t = time.time() -start
        
        self.matrix = matrix.T

        if verbose :
            print("Making matrix of size" , self.dim, "in", t, "sec.", "using core", core)

    def getMatrix(self):
        return self.matrix        
        
    def _isShift(self, fft=True):
        if self.dim % 2 == 0 and self.scaleinfo.domain[1][0]*self.scaleinfo.domain[1][1] < 0:
            isShift = True
        else:
            isShift = False
        if not fft:
            isShift = False
        return isShift
    
class SymetricUnitary(Unitary):
    def __init__(self, map, dim, tau):
        Unitary.__init__(self, map, dim, tau)

    def op0(self, x, isShift=False):
        if isShift:
            sub_x = numpy.fft.fftshift(x)
            func = self.map.ifunc0(sub_x)
            self.operator[0] = mpArray([ mpmath.exp(mpmath.mpc("0", -twopi*func[i]*self.tau/self.h*mpmath.mpf("0.5"))) for i in range(self.dim) ])
        else:
            func = self.map.ifunc0(x)
            self.operator[0] = mpArray([ mpmath.exp(mpmath.mpc("0", -twopi*func[i]*self.tau/self.h*mpmath.mpf("0.5"))) for i in range(self.dim) ])    

    def operate(self, invec=None, queue=None, verbose=False):
        if invec != None:
            if isinstance(invec, mpArray):
                self.setIn(invec)
            elif isinstance(invec, list) and self.dim==len(invec):
                self.setIn(mpArray(invec))

        qvec = self.operator[0]*self.stateIn
        pvec = mpArray(mpfft.fft(qvec, inverse=False, verbose=verbose))
        
        
        pvec = self.operator[1]*pvec #[self.operator[1].data[i]*pvec[i] for i in range(self.dim)]
        qvec = mpArray(mpfft.fft(pvec, inverse=True, verbose=verbose))
        
        qvec =self.operator[0]*qvec
        pvec = mpArray(mpfft.fft(qvec, inverse=False, verbose=verbose))
        
        self.stateOut = mpArray(mpfft.fft(pvec, inverse=True, verbose=verbose))
        
        if queue !=None: 
            queue.put(self.stateOut)

class NonUnitary(Unitary):
    pass


class QmapSystem(object):
    def __init__(self, map, dim, type, domain, **kwargs):
        self.map = map
        self.dim = dim
        self.type = type
        self._Setting(domain=domain,**kwargs)
        self.evals = None
        self.evecs = None
        self.matrix = None
        self._isSaveVector=True
        
        
    def _Setting(self, domain, **kwargs):
        tau = kwargs['tau'] if "tau" in kwargs else 1
        if self.type in [None, 'Unitary', 'uniraty','u','U'] :
            self.qmap  = Unitary(self.map, self.dim, tau)
        elif self.type in ["S","Symetric"]:
            self.qmap = SymetricUnitary(self.map, self.dim, tau)
        else:
            raise TypeError("type Error")
        self.qmap.setDomain(domain)
        self.scaleinfo = self.qmap.scaleinfo
    
    def setMatrix(self, **kwargs):
        haskey = self._hasKey('fft', 'core', 'verbose','left', 'solver','sort', **kwargs)
        fft = kwargs['fft'] if haskey[0] else True
        core = kwargs['core'] if haskey[1]  else 2
        verbose = kwargs['verbose'] if haskey[2] else False        
        self.qmap.setMatrix(fft=fft, core=core, verbose=verbose)
        self.matrix = self.qmap.getMatrix()

    def getMatrix(self):
        if self.matrix == None:
            self.setMtrix()
        return self.matrix

    def getEigen(self, **kwargs):
        haskey = self._hasKey('left', 'solver','sort', **kwargs)
        left = kwargs['left'] if haskey[0] else False
        solver = kwargs['solver'] if haskey[1] else 'qeispack'
        verbose = kwargs['verbose'] if haskey[2] else False                
        
        if self.matrix == None:
            self.setMatrix(**kwargs)
        if self.evals == None or self.evecs == None:
            self.evals, evecs = self.matrix.eigen(left, solver, verbose=verbose)
            self.evecs = [ State(self.scaleinfo, evec) for evec in evecs ]
        return (self.evals, self.evecs)

    def getEnergy(self,**kwargs):
        evals, evecs = self.getEigen(**kwargs)
        energy = mpArray([1.j*self.qmap.h*mpmath.log(val)/self.qmap.tau for val in evals])
        #energy = [x.real for x in ene]
        #decay = [x.imag for x in ene]
        return energy
    
    def sortEigen(self, val=None, order=True, index=None):
        if self.evals==None:  raise AttributeError("eigenvectors, eigenvalues is None")
        if index == None:
            index =[i[0] for i in sorted(enumerate(val), key=lambda x:x[1])]
        if order:
            self.evals = [self.evals[i] for i in index]
            self.evecs = [mpArray(self.evecs[i]) for i in index]
        else:
            self.evals = [self.evals[i] for i in index[::-1]]
            self.evecs = [mpArray(self.evecs[i]) for i in index[::-1]]
        
    
    def saveEigen(self, **kwargs):
        if self.evals ==None or self.evecs == None:
            self.getEigen(**kwargs)
            
        self.saveEvecs()
        self.saveEvals()
    
    def saveEvecs(self):
        if self.evecs == None: raise AttributeError("does not exists eigen-vectors")

        for i, evec in enumerate(self.evecs):
            f = open("eigen_qrep_%d.dat" % i, "w")
            x = self.qmap.scaleinfo.x[0]
            abs2 = evec.abs2()
            re = evec.real()
            im = evec.imag()
            self._annotate(f)
            f.write("# q-representation eigenvector, %d-th\n" % i )            
            f.write("# EIGENVALUE REAL %s\n" % self.evals[i].real)
            f.write("# EIGENVALUE IMAG %s\n" % self.evals[i].imag)
            f.write("# q, evec(q)*conj(evec(q)), Real[evec(q)], Imag[evec(q)]\n")            
            [f.write("%s %s %s %s \n" % (x[j], abs2[j], re[j], im[j])) for j in range(self.dim)]
            f.close()
        self._isSaveVector=True
    
    def saveEvals(self,title='eigen_vals.dat'):
        if self.evals == None: raise AttributeError("does not exists eigen-values")
        energy = self.getEnergy()
        with open(title, "w") as f:
            self._annotate(f)
            f.write("# index, Re[eval], Im[eval], energy, deday")
            for i, eval in enumerate(self.evals):
                ene = energy[i]
                f.write("%d %s %s %s %s\n" % (i, eval.real, eval.imag, ene.real,ene.imag))
    
    def _annotate(self, f):
        import datetime
        f.write("# DATE: %s\n" % datetime.datetime.now() )
        f.write("# DIM %d\n" % self.dim )
        f.write("# QMIN %s\n" % self.qmap.scaleinfo.domain[0][0])
        f.write("# QMAX %s\n" % self.qmap.scaleinfo.domain[0][1])
        f.write("# PMIN %s\n" % self.qmap.scaleinfo.domain[1][0])
        f.write("# PMAX %s\n" % self.qmap.scaleinfo.domain[1][1])
        
    def saveHsm(self, hsm_region=None,grid=[100,100],core=2,verbose=False):
        if not self._isSaveVector:
            self.saveEvecs()            
        if hsm_region == None:
            hsm_region = self.qmap.scaleinfo.domain
        from utility import parallel as p
        p.multi_hsm(hsm_region, grid,core=core,verbose=verbose)
        
            

    def _hasKey(self,*args,**kwargs):
        return [args[i] in kwargs for i in range(len(args))]
    