import subprocess
import matplotlib.pyplot as plt
import sys
from sphTools import *
from writeXS import *



def readInEigenvalue(filename):
    for line in reversed(list(open(filename))):
        if 'Combined k-effective' in line:
            return line.split()[3]

def plotThis(MC_modFlux,MC_fuelFlux):
    plt.plot(MC_modFlux[0],label='MC mod')
    plt.plot(MOC_modFlux[0],label='MOC mod')
    plt.plot(MC_fuelFlux[0],label='MC fuel')
    plt.plot(MOC_fuelFlux[0],label='MOC fuel')
    plt.legend(loc='best')
    plt.show()



if len(sys.argv) > 1:
    if sys.argv[1] == "full":
        subprocess.run(['python','../GettingXS/grid_3x3.py'])
        subprocess.run(['rm','geometry.xml','materials.xml','mgxs.h5',   \
                        'settings.xml','statepoint.100.h5','summary.h5', \
                        'tallies.out','tallies.xml'])
        print("\n\n")

#subprocess.run(['cp','XS.py','../MOC'])
#subprocess.run(['python','../MOC/grid.py','quiet','0'])
subprocess.run(['cp','XS.py','XS_original.py'])

from fluxMC import   *
MC_modFlux = [ MC_modFlux0, MC_modFlux1, MC_modFlux2, MC_modFlux3, MC_modFlux4, MC_modFlux5, MC_modFlux6, MC_modFlux7, MC_modFlux8 ]
MC_fuelFlux = [ MC_fuelFlux0, MC_fuelFlux1, MC_fuelFlux2, MC_fuelFlux3, MC_fuelFlux4, MC_fuelFlux5, MC_fuelFlux6, MC_fuelFlux7, MC_fuelFlux8 ]


nPins   = len(MC_modFlux)
nGroups = len(MC_modFlux[0])

MC_k = readInEigenvalue(filename='output')
print('MC SOLUTION',MC_k)

for i in range(9): MC_fuelFlux[i].reverse()
for i in range(9): MC_modFlux[i].reverse()
normalizeToUnity( nPins, MC_modFlux, MC_fuelFlux  )

#plotThis(MC_modFlux,MC_fuelFlux)

xsFileName = 'XS.py'


for i in range(10):

    # fluxMOC = MOC(XS)
    subprocess.run(['cp',xsFileName,'../MOC'])
    subprocess.run(['python','../MOC/grid.py','quiet',str(i)])
    exec("from fluxMOC"+str(i)+" import *")


    print('SPH Iter:',i,"   ",MOC_k)
    MOC_modFlux = [ MOC_modFlux0, MOC_modFlux1, MOC_modFlux2, MOC_modFlux3, MOC_modFlux4, MOC_modFlux5, MOC_modFlux6, MOC_modFlux7, MOC_modFlux8 ]
    MOC_fuelFlux = [ MOC_fuelFlux0, MOC_fuelFlux1, MOC_fuelFlux2, MOC_fuelFlux3, MOC_fuelFlux4, MOC_fuelFlux5, MOC_fuelFlux6, MOC_fuelFlux7, MOC_fuelFlux8 ]

    normalizeToUnity( nPins, MOC_modFlux, MOC_fuelFlux )

        
    # sph = SPH(fluxMOC,fluxMC)
    sph = calcSPH(nPins,nGroups,MOC_modFlux,MOC_fuelFlux,MC_modFlux,MC_fuelFlux)

    # update XS values
    exec("from XS"+str(i)+" import *")
    xsFileName = 'XS'+str(i+1)+'.py'
    writeXS( sph, nGroups, nPins, xsFileName,                                     \
             fuelTotal0, fuelAbsorption0, fuelNuFission0, fuelChi0, fuelScatter0, \
             modTotal0,  modAbsorption0,  modNuFission0,  modChi0,  modScatter0 )












subprocess.run(['cp','XS_original.py','XS.py'])







