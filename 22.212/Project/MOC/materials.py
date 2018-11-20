from fuelXS import *
from modXS  import *

class material:
    def __init__(self,SigmaT,SigmaF,SigmaS_matrix,chi,SigmaA):
        self.SigT = SigmaT
        self.SigF = SigmaF
        self.SigS_matrix = SigmaS_matrix
        self.chi  = chi
        self.SigA = SigmaA

def getScatteringIntoG(gNow,SMatrix,phi):
    return sum([SMatrix[g_from][gNow]*phi[g_from] for g_from in range(len(SMatrix))])

def getFissionIntoG(gNow,chi,SigF,phi):
    return sum([chi[gNow]*SigF[g_from]*phi[g_from] for g_from in range(len(SigF))])








modV = [material(modTotal0,modNuFission,modScatter0,modChi,absorMod0),
        material(modTotal1,modNuFission,modScatter1,modChi,absorMod1),
        material(modTotal2,modNuFission,modScatter2,modChi,absorMod2),
        material(modTotal3,modNuFission,modScatter3,modChi,absorMod3),
        material(modTotal4,modNuFission,modScatter4,modChi,absorMod4),
        material(modTotal5,modNuFission,modScatter5,modChi,absorMod5),
        material(modTotal6,modNuFission,modScatter6,modChi,absorMod6),
        material(modTotal7,modNuFission,modScatter7,modChi,absorMod7),
        material(modTotal8,modNuFission,modScatter8,modChi,absorMod8)]


fuelV = [material(fuelTotal0,fuelNuFission0,fuelScatter0,fuelChi0,absorFuel0),
         material(fuelTotal1,fuelNuFission1,fuelScatter1,fuelChi1,absorFuel1),
         material(fuelTotal2,fuelNuFission2,fuelScatter2,fuelChi2,absorFuel2),
         material(fuelTotal3,fuelNuFission3,fuelScatter3,fuelChi3,absorFuel3),
         material(fuelTotal4,fuelNuFission4,fuelScatter4,fuelChi4,absorFuel4),
         material(fuelTotal5,fuelNuFission5,fuelScatter5,fuelChi5,absorFuel5),
         material(fuelTotal6,fuelNuFission6,fuelScatter6,fuelChi6,absorFuel6),
         material(fuelTotal7,fuelNuFission7,fuelScatter7,fuelChi7,absorFuel7),
         material(fuelTotal8,fuelNuFission8,fuelScatter8,fuelChi8,absorFuel8)]

















