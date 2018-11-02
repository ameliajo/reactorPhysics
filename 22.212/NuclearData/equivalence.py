import matplotlib.pyplot as plt
from ultrafine import *



energies = np.linspace(1000,1100,500)
f1 = plt.figure(0)
for E in energies:
    SigmaF = slbw(E,sig_pot,atoms_per_b_cm)
    P_F_to_M = get_P_F_to_M(SigmaF,d)
    SigmaF_p = sig_pot*atoms_per_b_cm
    SigmaEq = ( P_F_to_M*SigmaF ) / ( 1.0 - P_F_to_M )
    phi = ((SigmaF_p/E) + (SigmaEq/E)) / (SigmaF + SigmaEq)

    plt.plot(E,phi,'bo',markersize=1)

f1.show()
input()







