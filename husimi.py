import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft

K = 0.6
M = 100
hbar = 2*np.pi/M

x0 = 1.5 
p0 = 0.3
n0 = 50
m0 = 40
n = np.arange(M)
m = np.arange(M)
x = 2*np.pi*n/M
p = 2*np.pi*m/M

def FFT(psi_x):
    psi_p = fft(psi_x)
    if np.linalg.norm(psi_p) != 0:
        psi_p = psi_p/np.linalg.norm(psi_p) # normalization
    return psi_p
def IFFT(psi_p):
    psi_x = ifft(psi_p)
    if np.linalg.norm(psi_x) != 0:
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
    if np.linalg.norm(psi0) != 0:
        psi0 = psi0/np.linalg.norm(psi0) # normalization
    return psi0

def evolution_step(psi, K):
    V = np.exp(-1j/hbar*K*np.cos(x))
    P = np.exp(-1j/(2*hbar)*p**2)

    psi = psi*V
    psi= FFT(psi)
    psi = psi* P
    psi = IFFT(psi)
    if np.linalg.norm(psi) != 0:
        psi = psi/np.linalg.norm(psi) # normalization
    return psi

def Husimi(x0,p0,psi):
    psi0 = create_psi_G(x0,p0)
    Q[:,n0] = FFT(np.conj(psi0)*psi)

    return np.dot(np.conj(psi0),psi)

psi0 = create_psi_G(x0,p0)
psi = psi0
    
for t in range(40):
    print(t)
        
    Q = np.zeros((M,M), dtype=complex)
    p0 = 0
    for j in range(M):
        x0 = 2*np.pi*j/M
        psi0 = create_psi_G(x0,p0)
        Q[:,j] = FFT(np.conj(psi0)*psi)

    n = np.arange(M)
    X, Y = np.meshgrid(n,n) # some convenient grid
    plt.pcolor( X,Y, np.abs(Q)**2, cmap='hot')
    plt.colorbar() # color scale#psi = evolution_step(psi, K)
    plt.show()
    psi = evolution_step(psi, K)
    print("xdxdddd")
'''
plt.plot(x,np.abs(psi)**2)
plt.title("Time "+str(i).zfill(3))
plt.xlabel("X")
plt.ylabel("Probability")
plt.savefig("wave_packet" + str(i).zfill(3) + ".png")
plt.close()            
'''
            
            
            
            
            
            
            
            
            
            