import numpy as np
from numpy import (fft, random)
from scipy.integrate import quad

def eisenstein_hu(Theta27, Omega0, h, **kwargs):
    Gamma = Omega0 * h
    
    def L0(q):
        return np.log(2*np.exp(1) + 1.8*q)
    
    def C0(q):
        return 14.2 + 731 / (1 + 62.5 * q)
    
    def T0(k):
        q = k * Theta27**2 / Gamma
        return L0(q) / (L0(q) + C0(q) * q**2)
    
    return T0

def cdm(ns, T0, A=1.0, **kwargs):
    return lambda k: A * k**ns * T0(k)**2

def W_th(y):
    return 3 / y**3 * (np.sin(y) - y * np.cos(y))

def norm_integrant(P, R):
    return lambda k: P(k) / (2 * np.pi**2) * W_th(R * k)**2 * k**2

class Cosmology(dict):
    def __init__(self, **kwargs):
        for key, val in kwargs.items():
            self[key] = val
        self.normalise()
            
    def __getattr__(self, key):
        return self[key]
    
    def __setattr__(self, key, val):
        self[key] = val
        
    def P(self, k):
        return cdm(T0=self.transfer_function(**self), **self)(k)

    @property
    def shape(self):
        return (self.N,) * 3
    
    @property
    def V(self):
        return self.L**3
    
    @property
    def freq(self):
        N = self.N
        L = self.L
        return fft.fftfreq(N, L/(2 * np.pi * N))
    
    @property
    def k(self):
        return self.freq[np.indices(self.shape)]
    
    @property
    def k_abs(self):
        return np.sqrt((self.k**2).sum(axis=0))

    def normalise(self, from_zero=False):
        self.A = 1.0
        self.A = self.sigma8**2 / quad(
            norm_integrant(self.P, 8.0),
            [(2*np.pi)/self.L, 0][from_zero],
            np.inf)[0]
        
Planck2018 = Cosmology(
    transfer_function=eisenstein_hu,
    h=0.674,
    ns=0.965,
    Omega0=1.000,
    sigma8=0.811,
    Theta27=2.7255/2.7,
    N=128,
    L=32.0,
    A=1.0)