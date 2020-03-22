

import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft


M = 1000
hbar = 2*np.pi/M

x0 = 1
p0 = 0.01
n0 = 50
m0 = 40
n = np.arange(M)
m = np.arange(M)
x = 2*np.pi*n/M
p = 2*np.pi*m/M

def FFT(psi_x):
    psi_p = fft(psi_x)
    psi_p = psi_p/np.linalg.norm(psi_p) # normalization
    return psi_p
def IFFT(psi_p):
    psi_x = ifft(psi_p)
    psi_x = psi_x/np.linalg.norm(psi_x) # normalization
    return psi_x

def create_psi_G(x0,p0):
    n0 = np.floor(x0/hbar)
    m0 = np.floor(p0/hbar)
    psi0 = np.zeros(M, dtype=complex) # complex vector
    s = np.zeros(M)
    for l in range(-4,5):
        s += np.exp( -(np.pi/M)*(n - n0 + l*M)**2)
    psi0 = np.exp(2j*np.pi *n* m0/M)*s
    psi0 = psi0/np.linalg.norm(psi0) # normalization
    return psi0

psi0_G_x = create_psi_G(x0,p0)    
plt.plot(x,np.abs(psi0_G_x)**2)
plt.axvline(x=x0)
plt.show()

psi0_G_p = FFT(psi0_G_x)
    
plt.plot(p,np.abs(psi0_G_p)**2)
plt.axvline(x=p0)
plt.show()
##################################################################################

#   evolution

##################################################################################
K = 1.1
V = np.exp(-1j/hbar*K*np.cos(x))
P = np.exp(-1j/(2*hbar)*p**2)
def evolution_step(psi):
    psi = psi*V
    psi= FFT(psi)
    psi = psi* P
    psi = IFFT(psi)
    psi = psi/np.linalg.norm(psi)
    return psi

psi0 = create_psi_G(x0,p0)
psi = psi0
for i in range(5):
    if i%100 == 0:
        print(i)
    psi = evolution_step(psi)
plt.plot(x,np.abs(psi)**2)
plt.show()

psi0 = create_psi_G(x0,p0)
psi = psi0
for i in range(100):
    print(i)
    plt.plot(x,np.abs(psi)**2)
    plt.title("Time "+str(i).zfill(3))
    plt.xlabel("X")
    plt.ylabel("Probability")
    plt.savefig("wave_packet" + str(i).zfill(3) + ".png")
    plt.close()            
    psi = evolution_step(psi)

            
            
            
            
            
            
            
            
            
            