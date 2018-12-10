from math import *

class Nuclide:
    def __init__(self,N,r,openMC_vals):
        self.N   = N
        self.pot = 4.0*pi*r*r
        self.sigT = []
        self.sigF = []
        self.sigA = []
        self.nuBar = []
        self.importFromOpenMC(openMC_vals)

    def addXS(self,xsVals):
        self.sigT.append(xsVals[0])
        self.sigF.append(xsVals[1])
        self.sigA.append(xsVals[2])
        self.nuBar.append(xsVals[3])

    def convertToMacro(self):
        self.SigT = [sig*1E-24*self.N  for sig in self.sigT]
        self.SigF = [sig*1E-24*self.N  for sig in self.sigF]
        self.SigA = [sig*1E-24*self.N  for sig in self.sigA]

    def importFromOpenMC(self,openMC_vals):
        self.openMC_nuSigF = openMC_vals[0]
        self.openMC_SigT   = openMC_vals[1]
        self.openMC_Chi    = openMC_vals[2]




        
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





