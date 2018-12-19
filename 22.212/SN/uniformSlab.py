import matplotlib.pyplot as plt
import numpy as np
from sweep import *

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




##############################################################################
# Define quadrature 
##############################################################################
from numpy.polynomial.legendre import leggauss
class GaussLegendre:
    def __init__(self, N):
        self.N = N
        self.mu, self.wgt = leggauss(N)
	
def runSN(numCells,width,SigT,SigS,src,methods,uniform):
    dx = width/numCells

    cells = [cell(i,SigT,SigS,src,dx,uniform) for i in range(numCells)]

    for method in methods:
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
                print(counter)
                plt.plot(x,phi,label=method+" "+str(numCells))
                converged = True
            oldPhi = phi[:]
            if counter > 1e3: 
                plt.plot(x,phi,label=method+" "+str(numCells))
                print(numCells,"is not stable")




S = GaussLegendre(2)
S = GaussLegendre(4)

numCells = 200
width = 20.0
SigT = 1.0
SigS = 0.1 # 0.1, 0.5, 0.99
src  = 1.0



methods = ['step','diamond','stepCharacteristic','linearDiscontinuous']

#phiGood = runSN(numCells,width,SigT,SigS,src,methods,True)
phiGood = runSN(numCells,width,SigT,SigS,src,methods,False)

plt.xlabel('distance (cm)')
plt.ylabel('avg. scalar flux')
plt.legend(loc='best')
plt.title('S'+str(S.N)+' method for '+str(width)+'cm slab, spatial convergence for scatter = '+str(SigS))
plt.show()



def getConvergence(numCells,method,SigT,SigS,src,width):
    nGood = 10000
    phiGood = runSN(nGood,width,SigT,SigS,src,[method])
    sumErrors = []
    for nCell in numCells:
        error = []
        sumError = 0.0
        phi = runSN(nCell,width,SigT,SigS,src,[method])
        for i in range(nCell):
            phi_Val = phi[i]
            goodVal = phiGood[int(i*nGood/nCell)]
            sumError += abs(goodVal-phi_Val)/goodVal
        
        sumErrors.append(sumError)
    plt.plot(numCells,sumErrors,label=method)
    return sumErrors
    
"""
numCells = [90,100,200,400,500,1000,2000,5000]
numCells = [90,100,200,400,700,1000]

sumErrors = getConvergence(numCells,'step',SigT,SigS,src,width)
sumErrors = getConvergence(numCells,'diamond',SigT,SigS,src,width)
sumErrors = getConvergence(numCells,'stepCharacteristic',SigT,SigS,src,width)
sumErrors = getConvergence(numCells,'linearDiscontinuous',SigT,SigS,src,width)
plt.xlabel('# Cells used')
plt.ylabel('% Error Relative to very fine run')
plt.legend(loc='best')

plt.title('S'+str(S.N)+' method for '+str(width)+'cm slab, spatial convergence for scatter = '+str(SigS))
    
plt.show()
"""




"""
plt.title('S'+str(S.N)+' method for '+str(width)+'cm slab, spatial convergence for '+methods[0])
#plt.title('S'+str(S.N)+' method for '+str(width)+'cm slab, using '+str(numCells)+' cells')
plt.xlabel('distance (cm)')
plt.ylabel('avg. scalar flux')
plt.legend(loc='best')
plt.show()



"""









