import subprocess
import matplotlib.pyplot as plt
import sys
from sphTools import *
from writeXS import *


if len(sys.argv) > 1:
    if sys.argv[1] == "full":
        subprocess.run(['python','../GettingXS/grid_3x3.py'])
        subprocess.run(['rm','geometry.xml','materials.xml','mgxs.h5',   \
                        'settings.xml','statepoint.100.h5','summary.h5', \
                        'tallies.out','tallies.xml'])
        print("\n\n")

subprocess.run(['cp','XS.py','XS_original.py'])

from fluxMC import   *

nPins, nGroups = len(MC_modFlux), len(MC_modFlux[0])

MC_k = readInEigenvalue(filename='output')
print('MC SOLUTION',MC_k)

for i in range(9): MC_fuelFlux[i].reverse()
for i in range(9): MC_modFlux[i].reverse()
normalizeToUnity( nPins, MC_modFlux, MC_fuelFlux  )


xsFile = 'XS'

# fluxMOC = MOC(XS)
subprocess.run(['cp',xsFile+'.py','../MOC/XS.py'])
subprocess.run(['python','../MOC/grid.py','quiet',str(i)])
exec("from fluxMOC"+str(i)+" import *")

print('MOC init    ',MOC_k)
MOC_k_vals = [MOC_k]


for i in range(500):

    normalizeToUnity( nPins, MOC_modFlux, MOC_fuelFlux )

    #diff = [abs(MC_modFlux[0][i]-MOC_modFlux[0][i]) for i in range(nGroups)]
    #plt.plot(diff,label='MOC mod'+str(i))
    
    #diff = [abs(MC_fuelFlux[0][i]-MOC_fuelFlux[0][i]) for i in range(nGroups)]
    #plt.plot(diff,label='MOC fuel'+str(i))
 

    # sph = SPH(fluxMOC,fluxMC)
    sph = calcSPH(nPins,nGroups,MOC_modFlux,MOC_fuelFlux,MC_modFlux,MC_fuelFlux)

    # update XS values
    exec("from "+xsFile+" import *")

    subprocess.run(['cp',xsFile+".py",'../MOC'])
    xsFile = 'XS'+str(i)
    writeXS( sph, nGroups, nPins, xsFile+".py",                                   \
             fuelTotal0, fuelAbsorption0, fuelNuFission0, fuelChi0, fuelScatter0, \
             modTotal0,  modAbsorption0,  modNuFission0,  modChi0,  modScatter0 )

    # fluxMOC = MOC(XS)
    subprocess.run(['cp',xsFile+'.py','../MOC/XS.py'])
    subprocess.run(['python','../MOC/grid.py','quiet',str(i)])
    exec("from fluxMOC"+str(i)+" import *")

    print('SPH Iter:',i,"   ",MOC_k)
    MOC_k_vals.append(MOC_k)


MC_k_vals = [MC_k]*len(MOC_k_vals)
plt.plot(MC_k_vals,label='MC k')
plt.plot(MOC_k_vals,label='MOC k')

plt.legend(loc='best')
plt.show()







subprocess.run(['cp','XS_original.py','XS.py'])







