from math import pi
import scipy
from scipy import special
import matplotlib.pyplot as plt
import numpy as np

sig_pot = 11.4
atoms_per_b_cm = 0.022
d = 0.4

def slbw(E,sig_pot,atoms_per_b_cm):
    T = 300.0
    Gamma = 0.095
    Gamma_n = 0.023
    Gamma_gamma = Gamma - Gamma_n
    E0 = 1050.0

    A = 238.0
    k = 8.6173303e-5

    r = (2603911.0/E0) * (A+1)/A
    x = 2.0 * (E-E0)/Gamma
    xsi = Gamma * ( A / ( 4*k*T*E0 ) )**0.5
    psi = (0.5 * (pi)**0.5 * xsi * scipy.special.wofz((x+1.0j)*xsi*0.5)).real

    sigma = (Gamma_n*Gamma_gamma)/(Gamma*Gamma) * (E0/E)**0.5 * r * psi
    Sigma = sigma*atoms_per_b_cm
    return Sigma

def get_P_F_to_M(SigmaF,d):
    E3 = scipy.special.expn(3,SigmaF*d)
    P_F_to_M = 1.0/(2.0*SigmaF*d) * (1.0 - 2.0*E3)
    return P_F_to_M


def transportCalc(E,sig_pot,atoms_per_b_cm,d):
    SigmaF_p = sig_pot*atoms_per_b_cm
    SigmaF = slbw(E,sig_pot,atoms_per_b_cm)
    P_F_to_M = get_P_F_to_M(SigmaF,d)
    phiF = ( (1.0-P_F_to_M)*(SigmaF_p/E)+P_F_to_M*(SigmaF/E) ) / SigmaF
    return phiF
 

def runUltraFine(sig_pot,atoms_per_b_cm,d):
    energies = np.linspace(1000,1100,500)
    #f1 = plt.figure(0)
    for E in energies:
        phiF = transportCalc(E,sig_pot,atoms_per_b_cm,d)
        SigmaF = slbw(E,sig_pot,atoms_per_b_cm)
        plt.plot(E,phiF,'ro',markersize=1)
        #plt.plot(E,SigmaF,'bo',markersize=1)
        #plt.plot(E,0.001*SigmaF,'bo',markersize=1)
    #f1.show()
    #input()











