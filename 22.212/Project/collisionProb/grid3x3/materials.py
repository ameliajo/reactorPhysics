import numpy as np
from XS import *


class Material:
    def __init__(self,SigT,nuSigF,SigS_matrix,chi,name,g):
        self.SigT = SigT
        self.nuSigF = nuSigF
        self.SigS = SigS_matrix
        self.chi = chi
        self.name = name

mMat  = Material(modTotal0,modNuFission0,modScatter0,modChi0,'mod',0)
fMat = Material(fuelTotal0,fuelNuFission0,fuelScatter0,fuelChi0,'fuel',0)



