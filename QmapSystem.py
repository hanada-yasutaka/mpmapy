from HilbertSpace import *
import mpfft
import numpy
import multiprocessing
twopi = mpmath.mpf("2")*mpmath.pi

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
            raise ValueError('Error: qmin >= qmax')
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

class State(Vector):
    def __init__(self, dim, scaleinfo, data=None):
        if not isinstance(scaleinfo, ScaleInfo):
            raise TypeError("scaleinfo must an instance with ScaleInfo")
        if data == None:
            Vector.__init__(self,input=dim)
        elif dim != len(data):
            raise ValueError("data length dose not much dim")
        else:
            Vector.__init__(self, input=data)
        
        self.scaleinfo = scaleinfo
        
    def set_qconst(self, q_c):
        q = self.scaleinfo.x[0]
        dx = mpmath.fabs(q[1] - q[0])
        qmin = self.scaleinfo.domain[0][0]
        qmax = self.scaleinfo.domain[0][1]
        index = int(round((q_c - qmin)/dx))
        if index > self.dim - 1 or index < 0:
            raise ValueError('q_c set in (%s,%s) with %s' %\
             (mpmath.nstr(qmin), mpmath.nstr(qmax), mpmath.nstr(dx)) )
        self.data[index] = mpmath.mpc(1, 0)
        
    def set_pconst(self, p_c):
        pass
    
    def set_cs(self,p_c):
        pass
    
    def set_linear(self):
        pass
        
    def qrep(self):
        return self.data
    
    def prep(self):
        pass
    
    def hsmrep(self, col, row, hsm_region=None):
        if hsm_region==None:
            self._hsm_region = self.scaleinfo.domain
        else:
            self._hsm_region = hsm_region 
        from utility import ctypes_wrapper
        cw = ctypes_wrapper.call_hsm_rep()
        X,Y, hsm_imag = cw.husimi_rep(self.data, self.dim, self.scaleinfo.domain, self._hsm_region, [row,col])        
        return X, Y, hsm_imag
    
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
        self.stateIn = Vector(dim)
        self.stateOut = Vector(dim)
        
    def setIn(self, state):
        if len(state) != self.dim and not isinstance(state, Vector):
            raise ValueError("expected State or Vector (dim=%d) instance" % self.dim)
        self.stateIn = state


    def getIn(self):
        return State(self.dim, self.scaleinfo, self.stateIn)
    
    def getOut(self):
        return State(self.dim, self.scaleinfo, self.stateOut)

    def swap(self):
        tmp = self.stateIn.data #[self.stateIn.data[i] for i in range(self.dim)]
        self.stateIn  = self.stateOut
        self.stateOut = Vector(tmp)

    def pull(self):
        self.stateIn = Vector(self.stateOut)

    def operate(self):
        pass
    
    def absorbed_operate(self):
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
        self.evals = None
        self.evecs = None
        
    def setDomain(self, r):
        self.scaleinfo.setDomain(r)
        self.h = self.scaleinfo.getPlanck()
        
    def op0(self, x, isShift=False):
        if isShift:
            sub_x = numpy.fft.fftshift(x)
            func = self.map.ifunc0(sub_x)
            self.operator[0] = Vector([ mpmath.exp(mpmath.mpc("0", -twopi*func[i]*self.tau/self.h)) for i in range(self.dim) ])
        else:
            func = self.map.ifunc0(x)
            self.operator[0] = Vector([ mpmath.exp(mpmath.mpc("0", -twopi*func[i]*self.tau/self.h)) for i in range(self.dim) ])

    def op1(self, x, isShift=True):
        if isShift:
            sub_x = numpy.fft.fftshift(x)
            func = self.map.ifunc1(sub_x)
            self.operator[1] = Vector([ mpmath.exp(mpmath.mpc("0", -twopi*func[i]*self.tau/self.h)) for i in range(self.dim) ])
        else:
            func = self.map.ifunc1(x)
            self.operator[1] = Vector([ mpmath.exp(mpmath.mpc("0", -twopi*func[i]*self.tau/self.h)) for i in range(self.dim) ])
        
    def operate(self, invec=None, queue=None, verbose=False):
        if invec != None:
            if isinstance(invec, Vector):
                self.setIn(invec)
            elif isinstance(invec, list) and self.dim==len(invec):
                self.setIn(Vector(invec))

        qvec = self.operator[0]*self.stateIn
        pvec = Vector(mpfft.fft(qvec.data, inverse=False, verbose=verbose))
        pvec = self.operator[1]*pvec #[self.operator[1].data[i]*pvec[i] for i in range(self.dim)]
        self.stateOut = Vector(mpfft.fft(pvec.data, inverse=True, verbose=verbose))
        if queue !=None: 
            queue.put(self.stateOut.data)
    

    def setOperator(self, fft=True):
        self.operator = [Vector(self.dim) for i in range(2)] 
        
        self.op0(self.scaleinfo.x[0], False)
        self.op1(self.scaleinfo.x[1], self._isShift())        

    def setMatrix(self, fft=True, core=2, verbose=False):
        from utility import parallel         
        self.setOperator(fft)
        self.matrix = Matrix(self.dim)

        orth_basis = [Vector(self.dim) for i in range(self.dim)]
        
        for i in range(self.dim):
            orth_basis[i].set(i, mpmath.mpc("1","0"))

        start = time.time()                            
        if core > 1:        
            data = [orth_basis[i].val() for i in range(self.dim)]            
            multi = parallel.MultiProcessing(self.operate, data, core=core)
            multi.run()
            res = multi.getResults()
            for i in range(self.dim):
                self.matrix.data[i,:] = mpmath.matrix(res[i]).T
        elif core == 1:            
            for i in range(self.dim):
                self.setIn(orth_basis[i])
                self.operate()
                self.matrix.data[i,:] = mpmath.matrix(self.stateOut.data).T
        t = time.time() -start
        
        self.matrix.data = self.matrix.data.T

        if verbose :
            print("Making matrix of size" , self.dim, "in", t, "sec.", "using core", core)

    def getMatrix(self):
        return self.matrix        
        
    def eigen(self, left=False, solver='qeispack', verbose=False):
        try:
            self.matrix.eigen(left=left, solver=solver, verbose=verbose)
        except AttributeError:
            self.setMatrix()
            self.matrix.eigen(left=left, solver=solver, verbose=verbose)
        self.evals, self.evecs = self.matrix.getEigen()
        self.evecs = [State(self.dim, self.scaleinfo, evec) for evec in self.evecs]

            
    def saveEigen(self,verbose=False):
        if self.evals ==None or self.evecs == None:
            self.getEigen(verbose=verbose)
        self._saveEvecs()
    
    def _saveEvals(self,title='eigen_vals.dat'):
        pass
    
    def _saveEvecs(self):
        for i, evec in enumerate(self.evecs):
            f = open("eigen_qrep_%d.dat" % i, "w")
            x = self.scaleinfo.x[0]
            abs2 = evec.abs2()
            self._annotate(f)
            f.write("# q-representation eigenvector, %d-th\n" % i )            
            f.write("# EIGENVALUE REAL %s\n" % self.evals.data[i].real)
            f.write("# EIGENVALUE IMAG %s\n" % self.evals.data[i].imag)
            f.write("# q, evec(q)*conj(evec(q)), Real[evec(q)], Imag[evec(q)]\n")            
            [f.write("%s %s %s %s \n" % (x[j], abs2[j], evec.data[j].real, evec.data[j].imag)) for j in range(self.dim)]
            f.close()

    def _annotate(self, f):
        import datetime
        f.write("# DATE: %s\n" % datetime.datetime.now() )
        f.write("# DIM %d\n" % self.dim )
        f.write("# QMIN %s\n" % self.scaleinfo.domain[0][0])
        f.write("# QMAX %s\n" % self.scaleinfo.domain[0][1])
        f.write("# PMIN %s\n" % self.scaleinfo.domain[1][0])
        f.write("# PMAX %s\n" % self.scaleinfo.domain[1][1])
        

    def _saveEigenHsm(self):
        pass       
            

    def _isShift(self, fft=True):
        if self.dim % 2 == 0 and self.scaleinfo.domain[1][0]*self.scaleinfo.domain[1][1] < 0:
            isShift = True
        else:
            isShift = False
        if not fft:
            isShift = False
        return isShift


class QmapSystem(object):
    def __init__(self, map, dim, type, domain, **kwargs):
        self.map = map
        self.dim = dim
        self.domain = domain
        self.type = type
        self._Setting()
        
    def _Setting(self):
        if self.type in [None, 'Unitary', 'uniraty','u','U'] :
            self.qmap  = Unitary(self.map, self.dim)
            self.qmap.setDomain(self.domain)
    
    def getEigen(self, **kwargs): #,fft=True, core=2, verbose=False):
        haskey = self._hasKey('fft', 'core', 'verbose','left', 'solver',**kwargs)
        fft = kwargs['fft'] if haskey[0] else True
        core = kwargs['core'] if haskey[1]  else 2
        verbose = kwargs['verbose'] if haskey[2] else False
        left = kwargs['left'] if haskey[3] else False
        solver = kwargs['left'] if haskey[3] else 'qeispack'        
        
        if self.qmap.evals == None or self.qmap.evecs == None:
            self.qmap.setMatrix(fft=fft, core=core, verbose=verbose)
            self.qmap.eigen(left, solver, verbose=verbose)
        return (self.qmap.evals, self.qmap.evecs)

    def getEnegry(self, **kwargs):
        eval, evecs = self.getEigen(**kwargs)
        qval = Vector([1.j*self.qmap.h*val.data/self.qmap.tau for val in eval])
        return qval
    
    def sorted(self):
        pass
    
    def saveEigen(self, **kwargs):
        haskey = self._hasKey('fft', 'core', 'verbose',**kwargs)
        fft = kwargs['fft'] if haskey[0] else True
        core = kwargs['core'] if haskey[1] else 2
        verbose = kwargs['verbose'] if haskey[2] else False

        if self.qmap.evals ==None or self.qmap.evecs == None:
            self.qmap.setMatrix(fft=fft, core=core, verbose=verbose)
            self.qmap.eigen(left=False, solver='qeispack',verbose=verbose)
            
        self.qmap.saveEigen(verbose=verbose)
    
    def saveEvecs(self):
        if self.qmap.evecs == None:
            raise ValueError("does not exists eigen-vectors") 
        for i, evec in enumerate(self.qmap.evecs):
            f = open("eigen_qrep_%d.dat" % i, "w")
            x = self.qmap.scaleinfo.x[0]
            abs2 = evec.abs2()
            self._annotate(f)
            f.write("# q-representation eigenvector, %d-th\n" % i )            
            f.write("# EIGENVALUE REAL %s\n" % self.qmap.evals.data[i].real)
            f.write("# EIGENVALUE IMAG %s\n" % self.qmap.evals.data[i].imag)
            f.write("# q, evec(q)*conj(evec(q)), Real[evec(q)], Imag[evec(q)]\n")            
            [f.write("%s %s %s %s \n" % (x[j], abs2[j], evec.data[j].real, evec.data[j].imag)) for j in range(self.dim)]
            f.close()
        self._isSave=True

    def _annotate(self, f):
        import datetime
        f.write("# DATE: %s\n" % datetime.datetime.now() )
        f.write("# DIM %d\n" % self.dim )
        f.write("# QMIN %s\n" % self.qmap.scaleinfo.domain[0][0])
        f.write("# QMAX %s\n" % self.qmap.scaleinfo.domain[0][1])
        f.write("# PMIN %s\n" % self.qmap.scaleinfo.domain[1][0])
        f.write("# PMAX %s\n" % self.qmap.scaleinfo.domain[1][1])
        
    def saveHsm(self, hsm_region=None,grid=[100,100],core=2,verbose=False):
        try:
            self._isSave
        except AttributeError:
            self.saveEvecs()
        if hsm_region == None:
            hsm_region = self.qmap.scaleinfo.domain
        from utility import parallel as p
        p.multi_hsm(hsm_region, grid,core=core,verbose=verbose)
        #for i,evec in enumerate(self.qmap.evecs):
        #    title = 'eigen_hsm_%d.dat' % i
        #    evec.save_hsmrep(col,row,hsm_region,title=title)
        
            
    def getMatrix(self,**kwargs):
        haskey = self._hasKey('fft', 'core', 'verbose',**kwargs)
        fft = kwargs['fft'] if haskey[0] else True
        core = kwargs['core'] if haskey[1] in kwargs else 2
        verbose = kwargs['verbose'] if haskey[2] else False

        if self.qmap.matrix == None:
            self.qmap.setMatrix(fft=fft, core=core, verbose=verbose)        
        return self.qmap.matrix
    
    def _hasKey(self,*args,**kwargs):
        return [args[i] in kwargs for i in range(len(args))]
    

        
    
            
        