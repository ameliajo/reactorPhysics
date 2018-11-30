import subprocess
import matplotlib.pyplot as plt
import sys


if len(sys.argv) > 1:
    if sys.argv[1] == "full":
        subprocess.run(['python','../GettingXS/grid_3x3.py'])
        subprocess.run(['rm','geometry.xml','materials.xml','mgxs.h5',   \
                        'settings.xml','statepoint.100.h5','summary.h5', \
                        'tallies.out','tallies.xml'])

subprocess.run(['cp','XS.py','../MOC'])
subprocess.run(['python','../MOC/grid.py','quiet','1'])

print("\n\n")

from nuclearData1 import *
from writeXS      import *
from sphTools     import *

nPins   = len(MC_modFlux)
nGroups = len(MC_modFlux[0])

MC_k = readInEigenvalue(filename='output')
for i in range(9): MC_fuelFlux[i].reverse()
for i in range(9): MC_modFlux[i].reverse()

print('MC SOLUTION',MC_k)
print('MOC Initial',MOC_k)

normalizeToUnity( nPins, MC_modFlux, MC_fuelFlux  )
normalizeToUnity( nPins, MOC_modFlux,MOC_fuelFlux )


plt.plot(MC_modFlux[0],label='MC mod')
plt.plot(MOC_modFlux[0],label='MOC mod')
plt.plot(MC_fuelFlux[0],label='MC fuel')
plt.plot(MOC_fuelFlux[0],label='MOC fuel')
plt.legend(loc='best')
plt.show()

sph = calcSPH(nPins,nGroups,MOC_modFlux,MOC_fuelFlux,MC_modFlux,MC_fuelFlux)
xsFileName = 'XS2.py'

print("\n")
print(sph[0])
print(sph[1])
print("\n")

writeXS( sph, nGroups, nPins, xsFileName,                                \
         fuelTotal, fuelAbsorption, fuelNuFission, fuelChi, fuelScatter, \
         modTotal,  modAbsorption,  modNuFission,  modChi,  modScatter )


subprocess.run(['cp',xsFileName,'../MOC/XS.py'])
subprocess.run(['python','../MOC/grid.py','quiet','2'])
from nuclearData2 import *

print('SPH Iter. 1',MOC_k)








normalizeToUnity( nPins, MC_modFlux, MC_fuelFlux  )
normalizeToUnity( nPins, MOC_modFlux,MOC_fuelFlux )

sph = calcSPH(nPins,nGroups,MOC_modFlux,MOC_fuelFlux,MC_modFlux,MC_fuelFlux)
xsFileName = 'XS3.py'

writeXS( sph, nGroups, nPins, xsFileName,                                \
         fuelTotal, fuelAbsorption, fuelNuFission, fuelChi, fuelScatter, \
         modTotal,  modAbsorption,  modNuFission,  modChi,  modScatter )


subprocess.run(['cp',xsFileName,'../MOC/XS.py'])
subprocess.run(['python','../MOC/grid.py','quiet','3'])
from nuclearData3 import *

print('SPH Iter. 3',MOC_k)







"""
for i in range(9):
    f.write("fuelTotal"+str(i)+"  = "+str([float("%.8f"%flux) for flux in phi[i][0]])+"\n")
    f.write("modTotal"+str(i)+" = "+str([float("%.8f"%flux) for flux in phi[i][1]])+"\n")
    f.write("\n\n")


"""



"""
#    plt.plot(MOC_mod,'ro-',label='MOC moderator'+str(i))
#    plt.plot(MC_mod, 'rx-',label='MC  moderator'+str(i))
#    plt.plot(MOC_fuel,'bo-',label='MOC fuel'+str(i))
#    plt.plot(MC_fuel, 'bx-',label='MC  fuel'+str(i))
#plt.legend(loc='best')
"""













