from fuelXS import *
from modXS  import *

class material:
    def __init__(self,SigmaT,SigmaF,SigmaS_matrix,chi):
        self.SigT = SigmaT
        self.SigF = SigmaF
        self.SigS_matrix = SigmaS_matrix
        self.chi  = chi

def getScatteringIntoG(gNow,SMatrix,phi):
    return sum([SMatrix[g_from][gNow]*phi[g_from] for g_from in range(len(SMatrix))])

def getFissionIntoG(gNow,chi,SigF,phi):
    return sum([chi[gNow]*SigF[g_from]*phi[g_from] for g_from in range(len(SigF))])








modV = [material(modTotal0,modNuFission,modScatter0,modChi),
        material(modTotal1,modNuFission,modScatter1,modChi),
        material(modTotal2,modNuFission,modScatter2,modChi),
        material(modTotal3,modNuFission,modScatter3,modChi),
        material(modTotal4,modNuFission,modScatter4,modChi),
        material(modTotal5,modNuFission,modScatter5,modChi),
        material(modTotal6,modNuFission,modScatter6,modChi),
        material(modTotal7,modNuFission,modScatter7,modChi),
        material(modTotal8,modNuFission,modScatter8,modChi)]


fuelV = [material(fuelTotal0,fuelNuFission0,fuelScatter0,fuelChi0),
         material(fuelTotal1,fuelNuFission1,fuelScatter1,fuelChi1),
         material(fuelTotal2,fuelNuFission2,fuelScatter2,fuelChi2),
         material(fuelTotal3,fuelNuFission3,fuelScatter3,fuelChi3),
         material(fuelTotal4,fuelNuFission4,fuelScatter4,fuelChi4),
         material(fuelTotal5,fuelNuFission5,fuelScatter5,fuelChi5),
         material(fuelTotal6,fuelNuFission6,fuelScatter6,fuelChi6),
         material(fuelTotal7,fuelNuFission7,fuelScatter7,fuelChi7),
         material(fuelTotal8,fuelNuFission8,fuelScatter8,fuelChi8)]

















