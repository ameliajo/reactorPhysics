import numpy as np
from XS import *


class Material:
    def __init__(self,SigT,nuSigF,SigS_matrix,chi,name):
        self.SigT = np.array(SigT)
        self.nuSigF = np.array(nuSigF)
        self.SigS = np.array(SigS_matrix)
        self.chi = np.array(chi)
        self.name = name

modClass  = Material(modTotal0,modNuFission0,modScatter0,modChi0,'mod')
fuelClass = Material(fuelTotal0,fuelNuFission0,fuelScatter0,fuelChi0,'fuel')



