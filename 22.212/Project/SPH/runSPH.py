import subprocess
import matplotlib.pyplot as plt
import sys


if len(sys.argv) > 1:
    if sys.argv[1] == "full":
        subprocess.run(['python','../GettingXS/grid_3x3.py'])

print("\n\n")
subprocess.run(['cp','XS.py','../MOC'])
subprocess.run(['rm','geometry.xml','materials.xml','mgxs.h5',   \
                'settings.xml','statepoint.100.h5','summary.h5', \
                'tallies.out','tallies.xml'])
subprocess.run(['python','../MOC/grid.py','quiet','1'])


from nuclearData1 import *
from writeXS      import *
from sphTools     import *

nPins   = len(MC_modFlux)
nGroups = len(MC_modFlux[0])

MC_k = readInEigenvalue(filename='output')
print('MC SOLUTION',MC_k)
print('MOC Iter. 1',MOC_k)


normalizeToUnity( nPins, MC_modFlux, MC_fuelFlux  )
normalizeToUnity( nPins, MOC_modFlux,MOC_fuelFlux )


sph = calcSPH(nPins,nGroups,MOC_modFlux,MOC_fuelFlux,MC_modFlux,MC_fuelFlux)
xsFileName = 'XS2.py'


writeXS( sph, nGroups, nPins, xsFileName,                                \
         fuelTotal, fuelAbsorption, fuelNuFission, fuelChi, fuelScatter, \
         modTotal,  modAbsorption,  modNuFission,  modChi,  modScatter )


subprocess.run(['cp',xsFileName,'../MOC/XS.py'])
subprocess.run(['python','../MOC/grid.py','quiet','2'])
from nuclearData2 import *

print('MOC Iter. 2',MOC_k)

 



"""
normalizeToUnity( nPins, MC_modFlux, MC_fuelFlux  )
normalizeToUnity( nPins, MOC_modFlux,MOC_fuelFlux )

sph = calcSPH(nPins,nGroups,MOC_modFlux,MOC_fuelFlux,MC_modFlux,MC_fuelFlux)
print("---------------------- SPH Val ",sph[0][0][0])

writeXS( sph, nGroups, nPins, 'XS3.py',                                  \
         fuelTotal, fuelAbsorption, fuelNuFission, fuelChi, fuelScatter, \
         modTotal,  modAbsorption,  modNuFission,  modChi,  modScatter )

subprocess.run(['cp','XS3.py','../MOC/XS.py'])
subprocess.run(['python','../MOC/grid.py','quiet','3'])
from nuclearData3 import *
print('MOC Iter. 3',MOC_k)
sph = calcSPH(nPins,nGroups,MOC_modFlux,MOC_fuelFlux,MC_modFlux,MC_fuelFlux)
print("---------------------- SPH Val ",sph[0][0][0])




subprocess.run(['rm','-rf','__pycache__'])

"""


















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













