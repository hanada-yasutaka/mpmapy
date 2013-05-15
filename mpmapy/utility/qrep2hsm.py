import os        
import sys
from ctypes_wrapper import call_hsm_rep
import numpy as np
import datetime


def grep(file, key):
    for line in open(file):
        if key in line:
            return line.split()[2]

def qrep2hsm(file, hsm_region, grid):
    qmin=float(grep(file, "QMIN"))
    qmax=float(grep(file, "QMAX"))
    pmin=float(grep(file, "PMIN"))
    pmax=float(grep(file, "PMAX"))
    dim=int(float(grep(file, "DIM")))
    data = np.loadtxt(file).transpose()
    qvec = data[2] + 1.j*data[3]
    cw = call_hsm_rep()

    hsm_imag = cw.husimi_rep(qvec, dim, [[qmin,qmax],[pmin,pmax]], hsm_region, grid)
    x = np.linspace(hsm_region[0][0], hsm_region[0][1], grid[0],endpoint=False)
    y = np.linspace(hsm_region[1][0], hsm_region[1][1], grid[1],endpoint=False)
    X,Y = np.meshgrid(x,y)
    data = np.array([X,Y,hsm_imag])

    #print hsm_imag[0,:]
#    print XX

    with open(file.replace("qrep","hsm"), "w") as of:
        of = open(file.replace("qrep","hsm"),"w") 
        of.write("# DATE %s\n" % datetime.datetime.now())
        of.write("# DIM %d\n# QMIN %f\n# QMAX %f\n" % (dim, qmin, qmax))
        of.write("# PMIN %f\n# PMAX %f\n" % (pmin, pmax))
        of.write("# VQMIN %f\n# VQMAX %f\n# VPMIN %f\n# VPMAX %f\n" % (hsm_region[0][0], hsm_region[0][1], hsm_region[1][0],hsm_region[1][1]))
        of.write("# ROW %d\n# COL %d\n" % (grid[0],grid[1]))
        for i,slice_data in enumerate(data.transpose()):
            np.savetxt(of, slice_data)
            of.write("\n")