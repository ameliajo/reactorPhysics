#from fuelXS import *
#from modXS  import *
from XS import *

class material:
    def __init__(self,SigmaT,SigmaF,SigmaS_matrix,chi,SigmaA):
        self.SigT = SigmaT
        self.SigF = SigmaF
        self.SigS_matrix = SigmaS_matrix
        self.chi  = chi
        self.SigA = SigmaA

def getScatteringIntoG(gNow,SMatrix,phi):
    #return sum([SMatrix[gNow][g_from]*phi[g_from] for g_from in range(len(SMatrix))])
    return sum([SMatrix[g_from][gNow]*phi[g_from] for g_from in range(len(SMatrix))])

def getFissionIntoG(gNow,chi,SigF,phi):
    return sum([chi[gNow]*SigF[g_from]*phi[g_from] for g_from in range(len(SigF))])








modV = [material(modTotal0,modNuFission0,modScatter0,modChi0,modAbsorption0),
        material(modTotal1,modNuFission1,modScatter1,modChi1,modAbsorption1),
        material(modTotal2,modNuFission2,modScatter2,modChi2,modAbsorption2),
        material(modTotal3,modNuFission3,modScatter3,modChi3,modAbsorption3),
        material(modTotal4,modNuFission4,modScatter4,modChi4,modAbsorption4),
        material(modTotal5,modNuFission5,modScatter5,modChi5,modAbsorption5),
        material(modTotal6,modNuFission6,modScatter6,modChi6,modAbsorption6),
        material(modTotal7,modNuFission7,modScatter7,modChi7,modAbsorption7),
        material(modTotal8,modNuFission8,modScatter8,modChi8,modAbsorption8)]


fuelV = [material(fuelTotal0,fuelNuFission0,fuelScatter0,fuelChi0,fuelAbsorption0),
         material(fuelTotal1,fuelNuFission1,fuelScatter1,fuelChi1,fuelAbsorption1),
         material(fuelTotal2,fuelNuFission2,fuelScatter2,fuelChi2,fuelAbsorption2),
         material(fuelTotal3,fuelNuFission3,fuelScatter3,fuelChi3,fuelAbsorption3),
         material(fuelTotal4,fuelNuFission4,fuelScatter4,fuelChi4,fuelAbsorption4),
         material(fuelTotal5,fuelNuFission5,fuelScatter5,fuelChi5,fuelAbsorption5),
         material(fuelTotal6,fuelNuFission6,fuelScatter6,fuelChi6,fuelAbsorption6),
         material(fuelTotal7,fuelNuFission7,fuelScatter7,fuelChi7,fuelAbsorption7),
         material(fuelTotal8,fuelNuFission8,fuelScatter8,fuelChi8,fuelAbsorption8)]

















