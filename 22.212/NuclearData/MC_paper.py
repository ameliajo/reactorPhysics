from math import pi
import scipy
from scipy import special
import matplotlib.pyplot as plt
import numpy as np

sig_pot = 11.4
atoms_per_b_cm = 0.022
d = 0.4
background = 8  # 8 b background fuel cross section
Sigma_background = background * atoms_per_b_cm
b = 1.15

energies = np.linspace(1000,1100,500)
#energies = np.linspace(1040,1060,100)

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
 

def runUltraFine(sig_pot,atoms_per_b_cm,d,plotMark):
    phi = [transportCalc(E,sig_pot,atoms_per_b_cm,d) for E in energies]
    plt.plot(energies,phi,plotMark,markersize=5,label="Ultrafine")
    





def runWigner(sig_pot,atoms_per_b_cm,d,Sigma_background,plotMark):
    phi = []
    for E in energies:
        l_bar = 2.0*d       # l_bar = 4 * V / S = 4 * (d*w*h) / (2*w*h) = 2*d
        SigmaEq = 1.0/l_bar + Sigma_background
        SigmaF = slbw(E,sig_pot,atoms_per_b_cm)
        SigmaF_p = sig_pot*atoms_per_b_cm
        phi.append(((SigmaF_p/E) + (SigmaEq/E)) / (SigmaF + SigmaEq))
    plt.plot(energies,phi,plotMark,markersize=3,label="Wigner")

def runBellWigner(sig_pot,atoms_per_b_cm,d,b,Sigma_background,plotMark):
    phi = []
    for E in energies:
        l_bar = 2.0*d       # l_bar = 4 * V / S = 4 * (d*w*h) / (2*w*h) = 2*d
        SigmaF = slbw(E,sig_pot,atoms_per_b_cm)
        #P_F_to_M = 1.0 / (SigmaF * l_bar + 1.0)
        SigmaEq = b/l_bar + Sigma_background
        SigmaF_p = sig_pot*atoms_per_b_cm
        phi.append(((SigmaF_p/E) + (SigmaEq/E)) / (SigmaF + SigmaEq))
    plt.plot(energies,phi,plotMark,markersize=3,label="BellWigner")


def runRoman(sig_pot,atoms_per_b_cm,d,b,Sigma_background,plotMark):
    alpha1 = 1.4
    alpha2 = 5.4
    beta   = 1.1
    phi = []
    for E in energies:
        SigmaE = b/(2.0*d)
        SigmaF_p = sig_pot*atoms_per_b_cm
        SigmaF = slbw(E,sig_pot,atoms_per_b_cm)
        P_F_to_M = beta * ( ( alpha1 * SigmaE ) / (SigmaF + alpha1*SigmaE) ) + \
             (1.0-beta) * ( ( alpha2 * SigmaE ) / (SigmaF + alpha2*SigmaE) ) 
        SigmaEq = ( ( P_F_to_M*SigmaF ) / ( 1.0 - P_F_to_M ) ) * Sigma_background
        phi.append(((SigmaF_p/E) + (SigmaEq/E)) / (SigmaF + SigmaEq))
    plt.plot(energies,phi,plotMark,markersize=3,label="Roman")



def runSubgroup(sig_pot,background):

    sigma_n = [0.0008, 0.0263, 0.1815, 2.4206, 9.9430]
    wgt = [0.3644, 0.3662, 0.2020, 0.0603, 0.0071]

    sigma_g_num = 0.0
    sigma_g_den = 0.0
    for i in range(len(sigma_n)):
        phi_sigma_n = (background + sig_pot)/(background+sigma_n[i])
        sigma_g_num += wgt[i]*sigma_n[i]*phi_sigma_n
        sigma_g_den += wgt[i]*phi_sigma_n

    return (sigma_g_num/sigma_g_den)



def runCarlvikWithDancoff(b,d,Sigma_background):
    C_vec = np.linspace(0.1,1.0,100)

    for E in np.linspace(1000,1100,9):
        sigma_eq_vec = []
        for C in C_vec:
            A = (1.0-C)/C
            alpha1 = (5.0*A+6.0-(A*A+36.0*A+36)**0.5)/(2.0*A+2.0)
            alpha2 = (5.0*A+6.0+(A*A+36.0*A+36)**0.5)/(2.0*A+2.0)
            beta = ((4.0*A+6.0)/(A+1.0) - alpha1)/(alpha2 - alpha1)

            SigmaE = b/(2.0*d)
            SigmaF = slbw(E,sig_pot,atoms_per_b_cm)

            P_F_to_M = beta * ( ( alpha1 * SigmaE ) / (SigmaF + alpha1*SigmaE) ) + \
                 (1.0-beta) * ( ( alpha2 * SigmaE ) / (SigmaF + alpha2*SigmaE) ) 
            SigmaEq = ( ( P_F_to_M*SigmaF ) / ( 1.0 - P_F_to_M ) ) * Sigma_background
            sigma_eq_vec.append(SigmaEq)
        plt.plot(C_vec,sigma_eq_vec,markersize=3,label="E = "+str(E))
    plt.title("Carlvik as a Function of Dancoff")
    plt.ylabel("SigmaEq")
    plt.xlabel("Dancoff")
    plt.legend(loc="upper right")
    plt.show()


runUltraFine(sig_pot,atoms_per_b_cm,d,'r')
runWigner(sig_pot,atoms_per_b_cm,d,Sigma_background,'y')
runBellWigner(sig_pot,atoms_per_b_cm,d,b,Sigma_background,'g')
runRoman(sig_pot,atoms_per_b_cm,d,b,Sigma_background,'k')
print("Subgroup Result",runSubgroup(sig_pot,background))
plt.xlabel("Energy (eV)")
plt.ylabel("Flux")
plt.show()



runUltraFine(sig_pot,atoms_per_b_cm,d,'r')
runWigner(sig_pot,atoms_per_b_cm,d,Sigma_background,'y')
runBellWigner(sig_pot,atoms_per_b_cm,d,b,Sigma_background,'g')
runRoman(sig_pot,atoms_per_b_cm,d,b,Sigma_background,'k')

plt.legend(loc="lower left")

d = 1.0
runUltraFine(sig_pot,atoms_per_b_cm,d,'r-.')
runWigner(sig_pot,atoms_per_b_cm,d,Sigma_background,'y-.')
runBellWigner(sig_pot,atoms_per_b_cm,d,b,Sigma_background,'g-.')
runRoman(sig_pot,atoms_per_b_cm,d,b,Sigma_background,'k-.')



plt.xlabel("Energy (eV)")
plt.ylabel("Flux")
plt.show()



runCarlvikWithDancoff(b,d,Sigma_background)








