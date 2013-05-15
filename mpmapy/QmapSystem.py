# -*- coding:utf-8 -*-
from mpArray import *
import mpfft
#import time
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
        self.x = [None] * 2
        self.h = None
        
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
        self.x = [mpArray(mpmath.linspace(self.domain[i][0], self.domain[i][1], self.dim, endpoint=False)) for i in range(2)]
        self.h = self.getPlanck()
        
    def getDomain(self):
        return self.domain

    def getX(self):
        return self.x

    def getVolume(self):
        return (self.domain[0][1] - self.domain[0][0]) * (self.domain[1][1] - self.domain[1][0])
    
    def getPlanck(self):
        return self.getVolume()/mpmath.mpf(self.dim)

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
    
    def savetxt(self, title, qp='q', abs2=False, header=None):
        with open(title, 'w') as f:
            if header != None:
                [f.write("%s" % head ) for head in header]
            if qp == 'p': x = self.scaleinfo.x[1]
            else: x = self.scaleinfo.x[0]
            re = self.real
            im = self.imag            
            if abs2:
                abs2 = self.abs2()
                [f.write("%s %s %s %s \n" % (x[j], abs2[j], re[j], im[j])) for j in range(len(self))]
            else:
                [f.write("%s %s %s \n" % (x[j], re[j], im[j])) for j in range(len(self))]
    
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

    def _cs(self, q_c, p_c, x=None):
        if not isinstance(q_c, mpmath.mpf) and not isinstance(p_c, mpmath.mpf):
            raise ValueError("q_c, p_c must be mpmath.mpf")
        if x == None:
            x = self.scaleinfo.x[0]
        re = -(x - q_c)*(x - q_c)*mpmath.pi/(self.scaleinfo.h)
        im = (x - q_c)*p_c*twopi/self.scaleinfo.h
        res = mpArray([mpmath.exp(mpmath.mpc(re[i], im[i])) for i in range(len(x))])
        return res.normalize()
    
    # todo: like classmethod 

    def cs(self, q_c, p_c):
        qrange = self.scaleinfo.domain[0]
        d = qrange[1] - qrange[0]
        lqmin, lqmax = qrange[0] - d, qrange[1] + d
        long_q = mpArray.linspace(lqmin, lqmax, 3*self.dim, endpoint=False)
        
        coh_state = self._cs(q_c, p_c, long_q)

        vec = mpArray(self.dim) 
        m = len(coh_state)/self.dim
        coh_state = coh_state.reshape(m,self.dim)
        
        for i in range(m):
            vec += coh_state[i][::1]
        return State(self.scaleinfo, vec.normalize())
        
           
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
            raise TypeError("tau types must be mpmath.mpf or integer")
    
        self.tau = tau
        self.matrix = None

        
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
        
        qvec = self.operator[0]*qvec
        pvec = mpArray(mpfft.fft(qvec, inverse=False, verbose=verbose))
        
        self.stateOut = mpArray(mpfft.fft(pvec, inverse=True, verbose=verbose))
        
        if queue !=None: 
            queue.put(self.stateOut)

class NonUnitary(Unitary):
    def __init__(self, map, dim, tau):
        Qmap.__init__(self, dim)
        if not isinstance(tau, mpmath.mpf) and  not isinstance(tau, int) and not isinstance(tau, mpmath.mpc):
            raise TypeError("tau types must be mpmath.mpf or mpmath.mpc or integer")
        
        self.map = map    
        self.tau = tau
        self.matrix = None        
        self.absorber = [mpArray.ones(dim) for i in range(2)]
        self.__dummy = [True, True,True]
        self.gamma = []
        
    def setAbsorber(self, qp=None, atype=None, region=None, **kwargs):
        start = time.time()
        qptypes = (None, 'q', 'p')
        atypes = ('box', 'cs', 'hole','exp', 'tanh','opt')

        if not qp in qptypes:
            raise TypeError("excepted first arg; 'qp' in", qptypes)
        
        if atype == 'box':
            absorber = self.box(qp, region)
        elif atype == 'hole' or atype == 'cs':
            absorber = self.hole(region, **kwargs)
            if self.__dummy[2]:
                self._operate = self.operate
                self.operate = self.hole_operate
                self.__dummy[2] = False
            qp = 'cs'
        elif atype == 'exp':
            absorber = self.abfunc(self.exp, qp, region, **kwargs)
        elif atype == 'tanh':
            absorber = self.abfunc(self.tanh, qp, region, **kwargs)
        elif atype == 'opt':
            absorber = self.abfunc(self.opt, qp, region, **kwargs)            
        else:
            absorber = None
            self.absorber = [ mpArray.ones(self.dim) for i in range(2)]
            print("Worning: absorbert does not set")

        self.addAbsorber(qp=qp, data = absorber, **kwargs)
                
        if True in (numpy.abs(self.absorber[0]) > 1.0) or True in (numpy.abs(self.absorber[1]) > 1.0):
            raise ValueError("souce term fund: absorber > 1.0:")

    
    def addAbsorber(self, qp, data=None, **kwargs):
        if data != None and qp == 'q':
            if self.__dummy[0]:
                self.absorber[0] = mpArray(self.dim)
                self.__dummy[0] = False

            self.absorber[0] += data

            if True in (numpy.abs(self.absorber[0]) > 1.0):
                self.absorber[0] = self.absorber[0]/self.absorber[0].real.max()            

            x = State(self.scaleinfo, self.absorber[0])
            x.savetxt("absorber_q.dat")
        
        elif qp =='p':
            if self.__dummy[1]:
                self.absorber[1] = mpArray(self.dim)
                self.__dummy[1] = False

            self.absorber[1] += data

            if True in (numpy.abs(self.absorber[1]) > 1.0):
                self.absorber[1] = self.absorber[1]/self.absorber[1].real.max()            

            x = State(self.scaleinfo, self.absorber[1])
            x.savetxt("absorber_p.dat",qp="p")            
        
        elif qp == 'cs':
            index = len(self.absorber)
            self.absorber.append(mpArray(self.dim))
            self.absorber[index] = data
            
            x = State(self.scaleinfo, mpArray.ones(self.dim) - self.absorber[index])
            x.savetxt("absorber_hole_%d.dat" % (index - 2))            
        
    def box(self, qp, region):
        if region == None: raise ValueError("region == None")
        if qp == 'q':
            index = (self.scaleinfo.x[0] > region[0]) & (self.scaleinfo.x[0] < region[1])
        elif qp == 'p':
            index = (self.scaleinfo.x[1] > region[0]) & (self.scaleinfo.x[1] < region[1])
        else:
            raise TypeError("excepted first arg is 'q' or 'p'")
        return mpArray.ones(self.dim)*index
        
    def hole(self, region, **kwargs):
        if not "gamma" in kwargs:
            raise ValueError("kwargs, gamma")
        self.gamma.append(kwargs['gamma'])

        return State(self.scaleinfo).cs(region[0], region[1])

    def abfunc(self, func, qp, region, **kwargs):
        if qp == 'q':
            abfunc = func(self.scaleinfo.x[0], region[0], region[1], **kwargs)
        elif qp == 'p':
            abfunc = func(self.scaleinfo.x[1], region[0], region[1], **kwargs)
        if 'cut' in kwargs:
            cut = kwargs['cut']
        else:
            cut = mpmath.mpf("1e-32")
        index = abfunc > cut
        return abfunc*index
    
    def opt(self, x, x1, x2, **kwargs):
        if not 'beta' in kwargs: raise ValueError("excepted kwargs: 'beta'")
        beta = kwargs['beta']
        expo = self._expotanh(x, x1, x2, **kwargs)
        return mpArray([mpmath.exp(expo[i]/(twopi*self.scaleinfo.h)*beta) for i in range(len(x))])
    
    def exp(self, x, x1,x2,**kwargs):
        if not 'beta' in kwargs: raise ValueError("excepted kwargs: 'beta'")
        beta = kwargs['beta']
        expo = self._expotanh(x, x1, x2, **kwargs)
        return mpArray([mpmath.exp(expo[i]/(twopi)*beta) for i in range(len(x))])
    
    def _expotanh(self, x, x1, x2, **kwargs):
        if 'alpha' in kwargs:
            alpha=kwargs['alpha']
        else:
            alpha=mpmath.mpf("1e10")
        print alpha
        z1 = x - x1
        z2 = x - x2
        y1 = mpArray([ mpmath.mpf("0.5")*(mpmath.mpf("1") - mpmath.tanh(x*alpha)) for x in z1 ]) 
        y2 = mpArray([ mpmath.mpf("0.5")*(mpmath.mpf("1") + mpmath.tanh(x*alpha)) for x in z2 ])

        return -(z1*z1*y1 + z2*z2*y2)        
        #return numpy.exp((-(z1**2*y1 + z2**2*y2)/(twopi*beta*h)))        
    
    def tanh(self, x, x1, x2, **kwargs):
        if not 'beta' in kwargs: raise ValueError("excepted kwargs: 'beta'")
        beta = kwargs['beta']
        z1 = (x - x1)*beta
        z2 = (x - x2)*beta
        return mpArray( [ mpmath.mpf("0.5")*(mpmath.tanh(z1[i]) - mpmath.tanh(z2[i])) for i in range(len(x)) ])

    def operate(self, invec=None, queue=None, verbose=False):
        if invec != None:
            self.setIn(invec)
        ## tood: 一度だけしか呼びたくない
        if self._isShift():
            self.absorber[1] = self.absorber[1].fftshift()
        qvec = self.absorber[0]*self.operator[0]*self.stateIn
        pvec = mpArray(mpfft.fft(qvec, inverse=False, verbose=verbose))
        pvec = self.absorber[1]*self.operator[1]*pvec 
        self.stateOut = self.absorber[0]*mpArray(mpfft.fft(pvec, inverse=True, verbose=verbose))
        
        if queue !=None: 
            queue.put(self.stateOut)
    
    def hole_operate(self, invec=None, queue=None, verbose=False):
        self._operate(invec, None, False) # = self.operate
        for i in range(2,len(self.absorber)):
            self.stateOut -= self.absorber[i]*self.stateOut.dot(self.absorber[i].conj())*self.gamma[i-2]*mpmath.mpf("0.5")
        if queue !=None:
            queue.put(self.stateOut)

class SymetricNonUnitary(NonUnitary):
    def __init__(self, map, dim, tau):
        NonUnitary.__init__(self, map, dim, tau)

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
        ## tood: 一度だけしか呼びたくない
        if self._isShift():
            self.absorber[1] = self.absorber[1].fftshift()

        qvec = self.absorber[0]*self.operator[0]*self.stateIn
        pvec = mpArray(mpfft.fft(qvec, inverse=False, verbose=verbose))
        
        pvec = self.absorber[1]*self.operator[1]*pvec
        qvec = mpArray(mpfft.fft(pvec, inverse=True, verbose=verbose))
        
        qvec = self.absorber[0]*self.operator[0]*qvec
        pvec = mpArray(mpfft.fft(qvec, inverse=False, verbose=verbose))
        
        self.stateOut = self.absorber[0]*mpArray(mpfft.fft(pvec, inverse=True, verbose=verbose))
        
        if queue !=None: 
            queue.put(self.stateOut)

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
        if self.type in [None, 'unitary' ,'u','U'] :
            self.qmap  = Unitary(self.map, self.dim, tau)
        elif self.type in ["S","symetric"]:
            self.qmap = SymetricUnitary(self.map, self.dim, tau)
        elif self.type in ['NU','nonunitary','open']:
            self.qmap = NonUnitary(self.map, self.dim, tau)
            self.setAbsorber = self.qmap.setAbsorber
        elif self.type in ['SNU', 'symetricnonunitary','snu']:
            self.qmap = SymetricNonUnitary(self.map, self.dim, tau)
            self.setAbsorber = self.qmap.setAbsorber
        else:
            raise TypeError("type Error")
        self.qmap.setDomain(domain)
        self.scaleinfo = self.qmap.scaleinfo
    
    def getState(self):
        return State(self.scaleinfo)
    
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
            re = evec.real
            im = evec.imag
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
    