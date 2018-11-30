
from XS1      import *
from fluxMC   import *
from fluxMOC1 import *


MC_modFlux = [ MC_modFlux0, MC_modFlux1, MC_modFlux2, MC_modFlux3, MC_modFlux4, MC_modFlux5, MC_modFlux6, MC_modFlux7, MC_modFlux8 ]
MC_fuelFlux = [ MC_fuelFlux0, MC_fuelFlux1, MC_fuelFlux2, MC_fuelFlux3, MC_fuelFlux4, MC_fuelFlux5, MC_fuelFlux6, MC_fuelFlux7, MC_fuelFlux8 ]

MOC_modFlux = [ MOC_modFlux0, MOC_modFlux1, MOC_modFlux2, MOC_modFlux3, MOC_modFlux4, MOC_modFlux5, MOC_modFlux6, MOC_modFlux7, MOC_modFlux8 ]
MOC_fuelFlux = [ MOC_fuelFlux0, MOC_fuelFlux1, MOC_fuelFlux2, MOC_fuelFlux3, MOC_fuelFlux4, MOC_fuelFlux5, MOC_fuelFlux6, MOC_fuelFlux7, MOC_fuelFlux8 ]


def readInEigenvalue(filename):
    for line in reversed(list(open(filename))):
        if 'Combined k-effective' in line:
            return line.split()[3]



