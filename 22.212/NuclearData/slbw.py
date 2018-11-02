from math import pi
import scipy
from scipy import special
import matplotlib.pyplot as plt
import numpy as np


def slbw(E):
    T = 300.0
    Gamma = 0.095
    Gamma_n = 0.023
    Gamma_gamma = Gamma - Gamma_n
    E0 = 1050.0
    sig_pot = 11.4

    A = 238.0
    k = 8.6173303e-5

    sig_background = 8.0
    SigmaMod = 1.23

    r = (2603911.0/E0) * (A+1)/A
    xsi = Gamma * ( A / ( 4*k*T*E0 ) )**0.5
    x = 2.0 * (E-E0)/Gamma
    q = 2.0 * (r*sig_pot)**0.5

    psi = pi**0.5 * (xsi * scipy.special.wofz((x+1.0j)*xsi)).real
    chi = pi**0.5 * (xsi * scipy.special.wofz((x+1.0j)*xsi)).imag

    sigGamma = (Gamma_n/Gamma)*(Gamma_gamma/Gamma)*(E0/E)**0.5*(r*psi)
    sigN     = (Gamma_n/Gamma)*(Gamma_n/Gamma)*(r*psi + (Gamma/Gamma_n)*q*chi)+sig_pot
    return sigGamma,sigN

print(slbw(3.5))


energies = np.linspace(1000,1100,1000)
f1 = plt.figure(0)
for E in energies:
    sigGamma = slbw(E)[0]
    plt.plot(E,sigGamma,'ro',markersize=1)
    sigN = slbw(E)[1]
    plt.plot(E,sigN,'bo',markersize=1)
f1.show()
input()


















