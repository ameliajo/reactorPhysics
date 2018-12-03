import matplotlib.pyplot as plt
import numpy as np
plt.style.use('seaborn-dark')

import openmc
import openmc.mgxs as mgxs
import urllib.request


def getPointwiseSig(E):
    return total.xs['294K'](E)


def fluxEquivalentEq(E, C, fuelSigP, SigEq):
    fuelSig = getPointwiseSig(E)
    return (fuelSigP + SigEq) / (E * (fuelSig + SigEq))




# Download ENDF file
url = 'https://t2.lanl.gov/nis/data/data/ENDFB-VII.1-neutron/U/235'
filename, headers = urllib.request.urlretrieve(url, 'u235.endf')

# Load into memory
u235 = openmc.data.IncidentNeutron.from_endf(filename)
total = u235[1]



C = 0.1452
fuelSigP = 2.0
SigEq = 1.0/(3.1)

#phi = []
#for E in range(1,1000):
#    phi.append(fluxEquivalentEq(E,C,fuelSigP,SigEq))
#plt.plot(phi)
#plt.show()



openmc.plot_xs('U235', ['total'])

# Download ENDF file
url = 'https://t2.lanl.gov/nis/data/data/ENDFB-VII.1-neutron/U/235'
filename, headers = urllib.request.urlretrieve(url, 'u235.endf')

# Load into memory
u235 = openmc.data.IncidentNeutron.from_endf(filename)
total = u235[1]
print(u235)

u235Total = []

E_vec = np.linspace(1e-5,1e6,1e6)
for E in E_vec:
    u235Total.append(total.xs['0K'](E))

plt.loglog(E_vec,u235Total,'r')
plt.ylim(1e-1,1e5)


print("finished plotting")


plt.show()







