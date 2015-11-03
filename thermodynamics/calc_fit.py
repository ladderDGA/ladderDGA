#!/bin/python

import numpy as np
import sys
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

name = sys.argv[1]
print name

w = np.loadtxt('w.dat',unpack=True)
freq,ek1_re,ek1_im,ek2_re,ek2_im,ep1_re,ep1_im,ep2_re,ep2_im = np.loadtxt(name,unpack=True)

fitpoints = np.size(w)
iwmax = np.size(freq)

ek1 = 0.0
ek2 = 0.0
ep1 = 0.0
ep2 = 0.0

start = iwmax - np.size(w)

for i in range(fitpoints):
  ek1=ek1+ek1_re[iwmax-fitpoints+i]*w[i]
  ek2=ek2+ek2_re[iwmax-fitpoints+i]*w[i]
  ep1=ep1+ep1_re[iwmax-fitpoints+i]*w[i]
  ep2=ep2+ep2_re[iwmax-fitpoints+i]*w[i]

print ek1,ek2,ep1,ep2
