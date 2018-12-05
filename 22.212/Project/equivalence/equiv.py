import matplotlib.pyplot as plt
import numpy as np
plt.style.use('seaborn-dark')
import openmc
import urllib.request


def getPointwiseSig(E, total):
    return total.xs['0K'](E)



def calcSigEq(approx,lbar,b=1):
    return 1.0/lbar if approx == 'wigner' else b/lbar


def fluxEquivalentEq(E, C, fuelSigP, SigEq, total):
    fuelSig = getPointwiseSig(E, total)
    return (fuelSigP + SigEq) / (E * (fuelSig + SigEq))



url = 'https://t2.lanl.gov/nis/data/data/ENDFB-VII.1-neutron/U/235'
filename, headers = urllib.request.urlretrieve(url, 'u235.endf')
u235 = openmc.data.IncidentNeutron.from_endf(filename)
total = u235[1]



C = 0.1452
fuelSigP = 2.0
lbar = 3.1

phi = []
for E in np.linspace(1e-1,1e3,1e2):
    SigEq = calcSigEq('wigner',lbar)
    phi.append(fluxEquivalentEq(E,C,fuelSigP,SigEq,total))
plt.plot(phi)
plt.show()






