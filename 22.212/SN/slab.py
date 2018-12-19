import matplotlib.pyplot as plt
import numpy as np
from sweep import *
from numpy.polynomial.legendre import leggauss
import sys



numCells = 200
width = 25.0
SigT = 1.0
SigS = 0.1 # 0.1, 0.5, 0.99
src  = 1.0
N = 2

flux = True
error = False
uniform = False



###############################################################################
class GaussLegendre:
    def __init__(self, N):
        self.N = N
        self.mu, self.wgt = leggauss(N)

S = GaussLegendre(N)


###############################################################################
class cell:
    def __init__(self,ID,SigT,SigS,src,width,uniform):
        self.ID    = ID
        self.SigT  = SigT
        self.SigS  = SigS
        self.src   = src
        self.phi   = 1.0
        self.dx    = width
        self.L     = ID*width
        self.M     = (ID+0.5)*width
        self.R     = (ID+1)*width
        if not uniform:
            if self.M <= 10:
                self.src = 0.0
                self.SigS = 0.5
            if 10 < self.M <= 20:
                self.src = 1.0
                self.SigS = 0.9
            if 20 < self.M:
                self.src = 0.0
                self.SigS = 0.2
    def __str__(cell):
        return "Cell #"+str(cell.ID)+" and spans "+str(cell.L)+" - "+str(cell.R)





def runSN(numCells,width,SigT,SigS,src,method,uniform,plot):
    dx = width/numCells

    cells = [cell(i,SigT,SigS,src,dx,uniform) for i in range(numCells)]

    oldPhi = [1.0]*numCells
    converged = False
    counter = 0
    while not converged:
        counter += 1
        broom(S,cells,method,'vacuum')
        x   = [cell.M   for cell in cells]
        phi = [cell.phi for cell in cells]
        invMaxVal = 1.0/max(phi)
        phi = [x*invMaxVal for x in phi]

        diff = [oldPhi[i]-phi[i] for i in range(numCells)]
        if max(diff) < 1e-6: 
            print(method,'using',numCells,'cells took',counter,'iterations to converge')
            if plot: plt.plot(x,phi,label=method+" "+str(numCells))
            converged = True
        oldPhi = phi[:]
        if counter > 1e3: 
            print(method,'using',numCells,"cells is not stable")

    return phi

def getConvergence(method,SigT,SigS,src,width):
    nGood = 10000
    numCells = [100,200,400,700,1000,5000,8000]
    phiGood = runSN(nGood,width,SigT,SigS,src,method,False,False)
    sumErrors = []
    for nCell in numCells:
        error = []
        sumError = 0.0
        phi = runSN(nCell,width,SigT,SigS,src,method,False,False)
        for i in range(nCell):
            phi_Val = phi[i]
            goodVal = phiGood[int(i*nGood/nCell)]
            sumError += abs(goodVal-phi_Val)/goodVal
        
        sumErrors.append(sumError)
    plt.plot(numCells,sumErrors,label=method)
 


if flux:
    methods = ['step','diamond','stepCharacteristic','linearDiscontinuous']
    phiGood = runSN(numCells,width,SigT,SigS,src,methods[0],uniform,True)
    phiGood = runSN(numCells,width,SigT,SigS,src,methods[1],uniform,True)
    phiGood = runSN(numCells,width,SigT,SigS,src,methods[2],uniform,True)
    phiGood = runSN(numCells,width,SigT,SigS,src,methods[3],uniform,True)

    plt.xlabel('distance (cm)')
    plt.ylabel('avg. scalar flux')
    plt.legend(loc='best')
    if uniform:
        plt.title('S'+str(S.N)+' method for '+str(width)+'cm slab, spatial convergence for scatter = '+str(SigS))
    else:
        plt.title('S'+str(S.N)+' method for '+str(width)+'cm slab, heterogeneous problem, '+str(numCells)+' cells')
    plt.show()



   
if error:
    getConvergence('step',SigT,SigS,src,width)
    getConvergence('diamond',SigT,SigS,src,width)
    getConvergence('stepCharacteristic',SigT,SigS,src,width)
    getConvergence('linearDiscontinuous',SigT,SigS,src,width)
    plt.xlabel('# Cells used')
    plt.ylabel('% Error Relative to very fine run')
    plt.legend(loc='best')
    if uniform: plt.title('S'+str(S.N)+' method for '+str(width)+\
                  'cm slab, spatial convergence for scatter = '+str(SigS))
    else:       plt.title('S'+str(S.N)+' method for '+str(width)+\
                  'cm slab, spatial convergence for hetero slab') 
    plt.show()











