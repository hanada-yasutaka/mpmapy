#!/usr/bin/env python 

from __future__ import division
import mpmath
import time


def dft(x, inverse = False, verbose = False) :
    t = time.time()
    N = len(x)
    inv = mpmath.mpf("-1") if not inverse else mpmath.mpf("1")
    X =[mpmath.mpf("0")] * N
    for k in range(N) :
        for n in range(N) :
            X[k] += x[n] * mpmath.exp(inv * mpmath.mpc("0","2") * mpmath.pi * mpmath.mpf(k) * mpmath.mpf(n) / mpmath.mpf(N))            
        if inverse :
            X[k] /= N
    t = time.time() - t
    if verbose :
        print("Computed","an inverse" if inverse else "a","DFT of size",N,)
        print("in",t,"sec.")
    return X


def fft(x, inverse = False, verbose = False) :
    t = time.time()
    N = len(x)
    inv = mpmath.mpf("-1") if not inverse else mpmath.mpf("1")
    if N % 2 != 0 :
        return dft(x, inverse, verbose)
    x_e = x[::2]
    x_o  = x[1::2]
    X_e = fft(x_e, inverse, False)
    X_o  = fft(x_o, inverse, False)
    X = []
    M = N // 2
    for k in range(M) :
        X += [X_e[k] + X_o[k] * mpmath.exp(inv * mpmath.mpc("0","2") * mpmath.pi * mpmath.mpf(k) / mpmath.mpf(N))]        
    for k in range(M,N) :
        X += [X_e[k-M] - X_o[k-M] * mpmath.exp(inv * mpmath.mpc("0","2") * mpmath.pi * mpmath.mpf(k-M) / mpmath.mpf(N))]
    if inverse :
        X = [mpmath.fdiv(j, mpmath.mpf(2)) for j in X]
    t = time.time() - t
    if verbose :
        print("Computed","an inverse" if inverse else "a","CT FFT of size",N,)
        print("in",t,"sec.")
    return X
