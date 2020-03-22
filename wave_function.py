# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 16:18:20 2020

@author: Antek
"""

import numpy as np
import matplotlib.pyplot as plt

M = 100
hbar = 2*np.pi/M

x0 = 1
n0 = 3
m0 = 10
n = np.arange(M)
m = np.arange(M)
x = 2*np.pi*n/M

psi0 = np.zeros(M, dtype=complex) # complex vector
psi0 = np.exp( -(x-x0)**2/2/hbar ) # element-wise
# x -- vector
psi0 = psi0/np.linalg.norm(psi0) # normalization
plt.plot(x,np.abs(psi0)**2)
plt.title("Wave function, x0 = "+str(x0))
plt.xlabel("X")
plt.ylabel("Probability")
plt.savefig("Wave_func.png")
