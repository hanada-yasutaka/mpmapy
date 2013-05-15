from multiprocessing import Queue, Process

class MultiProcessing(object):
    def __init__(self, f, data, core=2, **kwargs):
        self.f = f
        self.core=core
        self.kwargs=kwargs
        self.data_num = len(data)
        self.res = [None] * self.data_num        
        self.data = data
        #self.rum = self.samefunc
    
    def run(self):

        iter_max = int(self.data_num/self.core)
        k = 0
        for j in range(iter_max):
            queues = [ Queue() for i in range(self.core) ]

            p = [Process(target=self.f, args=(self.data[i+k],queues[i])) for i in range(self.core)]
            [p[i].start() for i in range(self.core)]
            for i in range(self.core):
                self.res[i+k] = queues[i].get()
            [p[i].join() for i in range(self.core)]            
            k = k + self.core
        rest = (self.data_num % self.core)
        if rest !=0:
            queues = [ Queue() for i in range(rest) ]
#            data_list = [managers[i].list(self.data[i+k]) for i in range(rest)]
            p = [Process(target=self.f, args=(self.data[i+k],queues[i])) for i in range(rest)]
            [p[i].start() for i in range(rest)]
            for i in range(rest):
                index = i+k
                self.res[index] = queues[i].get()
            [p[i].join() for i in range(rest)]
                        
    
    def getResults(self):
        if None in self.res:
            raise ValueError("can't return")
        return self.res 

def multi_hsm(hsm_region, grid, core, verbose=False):
    import time
    import qrep2hsm
    import sort
    import glob    
    start = time.time()

    files = glob.glob("eigen_qrep*.dat")
    data_num = len(files)
    sort.natural(files)
    iter_max = int(data_num/core)
    k = 0
    for j in range(iter_max):
        if verbose:
            print([(files[i+k],"->",files[i+k].replace("qrep","hsm")) for i in range(core)])
        #title="eigen_hsm_%d.dat" % j
        p = [Process(target=qrep2hsm.qrep2hsm, args=(files[i+k],hsm_region,grid)) for i in range(core)]
        [p[i].start() for i in range(core)]
        [p[i].join() for i in range(core)]            
        k = k + core
    rest = (data_num % core)
    if rest !=0:
        if verbose:
            print([(files[i+k],"->",files[i+k].replace("qrep","hsm")) for i in range(rest) ])
        p = [Process(target=qrep2hsm.qrep2hsm, args=(files[i+k],hsm_region,grid)) for i in range(rest)]
        [p[i].start() for i in range(rest)]
        for i in range(rest):
            index = i+k
        [p[i].join() for i in range(rest)]
    t = time.time() - start
    if verbose:
        print("number of",data_num,"qrep data to hsm data in", t,"sec. using core", core)       