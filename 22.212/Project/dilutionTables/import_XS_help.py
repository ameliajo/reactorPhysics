from math import *

class Nuclide:
    def __init__(self,N,r,openMC_vals):
        self.N   = N
        self.pot = 4.0*pi*r*r
        self.importFromOpenMC(openMC_vals)

    def resetXS(self,nGroups):
        self.sigT = [None]*nGroups
        self.sigF = [None]*nGroups
        self.sigA = [None]*nGroups
        self.nuBar = [None]*nGroups

    def addXS(self,xsVals,g):
        self.sigT[g] = (xsVals[0])
        self.sigF[g] = (xsVals[1])
        self.sigA[g] = (xsVals[2])
        self.nuBar[g] = (xsVals[3])

    def convertToMacro(self,nGroups):
        self.SigT = [self.sigT[g]*1E-24*self.N if self.sigT[g] else self.openMC_SigT[g] for g in range(nGroups)]
        self.SigA = [self.sigA[g]*1E-24*self.N if self.sigA[g] else self.openMC_SigA[g] for g in range(nGroups)]
        self.nuSigF = [self.nuBar[g]*self.sigF[g]*1E-24*self.N if self.sigF[g] else self.openMC_nuSigF[g] for g in range(nGroups)]

    def importFromOpenMC(self,openMC_vals):
        self.openMC_nuSigF = openMC_vals[0]
        self.openMC_SigT   = openMC_vals[1]
        self.openMC_Chi    = openMC_vals[2]
        self.openMC_SigA   = openMC_vals[3]
        self.openMC_SigS   = openMC_vals[4]

        self.nuSigF = openMC_vals[0]
        self.SigT   = openMC_vals[1]
        self.Chi    = openMC_vals[2]
        self.SigA   = openMC_vals[3]
        self.SigS   = openMC_vals[4]






        
def interpolate(x1,y1,x2,y2,x):
    return y1 + ((y2-y1)/(x2-x1))*x - (y2*x1-y1*x1)/(x2-x1)

def getDataFromEqTable(boundingDilutions,nuclideData,backgroundXS):
    d1     = boundingDilutions[0][1]
    sigT1  = nuclideData.sigT[boundingDilutions[0][0]]
    sigF1  = nuclideData.sigF[boundingDilutions[0][0]]
    sigA1  = nuclideData.sigA[boundingDilutions[0][0]]
    nuBar1 = nuclideData.nuBar
    
    d2     = boundingDilutions[1][1]
    sigT2  = nuclideData.sigT[boundingDilutions[1][0]]
    sigF2  = nuclideData.sigF[boundingDilutions[1][0]]
    sigA2  = nuclideData.sigA[boundingDilutions[1][0]]
    nuBar2 = nuclideData.nuBar

    sigT = (interpolate(d1,sigT1,d2,sigT2,backgroundXS))
    sigF = (interpolate(d1,sigF1,d2,sigF2,backgroundXS))
    sigA = (interpolate(d1,sigA1,d2,sigA2,backgroundXS))
    nuBar = (interpolate(d1,nuBar1,d2,nuBar2,backgroundXS))

    return sigT,sigF,sigA,nuBar





