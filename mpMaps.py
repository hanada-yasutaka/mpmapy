import mpmath
twopi = mpmath.mpf("2")*mpmath.pi

class StandardMap(object):
    def __init__(self, k):
        if not isinstance(k,mpmath.mpf) and not isinstance(k,int):
            raise TypeError("parameter k is expected mpmath.mpf or int type")
        self.k = k
        
    def ifunc0(self, x):
        return [self.k*mpmath.cos( twopi * x[i] )/(twopi*twopi) for i in range(len(x))]
    
    def ifunc1(self, x):
        return [x[i]*x[i]/mpmath.mpf("2") for i in range(len(x))]
        
if __name__ == '__main__':
    mpmath.mp.dps=40
    mpmath.mp.pretty=True
    map = StandardMap(2)
    x = mpmath.linspace(0, 1, 10)
    print(x)
    print(map.ifunc0(x))
    print(map.ifunc1(x))
        