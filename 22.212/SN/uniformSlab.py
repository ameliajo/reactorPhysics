import matplotlib.pyplot as plt
import numpy as np
from sweep import *

class cell:
    def __init__(self,ID,SigT,SigS,src,width):
        self.ID    = ID
        self.SigT  = SigT
        self.SigS  = SigS
        self.src   = src
        self.phi   = 1.0
        self.dx    = width
        self.L     = ID*width
        self.M     = (ID+0.5)*width
        self.R     = (ID+1)*width
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
	
def runSN(numCells,width,SigT,SigS,src,methods):
    dx = width/numCells
    cells = []
    for i in range(numCells):
        cells.append(cell(i,SigT,SigS,src,dx))



    for method in methods:
        oldPhi = [1.0]*numCells
        converged = False
        counter = 0
        while not converged:
            counter += 1
            broom(S,cells,method)
            x   = [cell.M   for cell in cells]
            phi = [cell.phi for cell in cells]
            invMaxVal = 1.0/max(phi)
            phi = [x*invMaxVal for x in phi]
            #plt.plot(x,phi,label=method)

            diff = [oldPhi[i]-phi[i] for i in range(numCells)]
            if max(diff) < 1e-6: 
                #plt.plot(x,phi,label=method+" "+str(numCells))
                converged = True
                return phi
            oldPhi = phi[:]
            if counter > 1e3: 
                #plt.plot(x,phi,label=method+" "+str(numCells))
                print(numCells,"is not stable")
                return phi




S = GaussLegendre(2)
#S = GaussLegendre(4)

numCells = 10000
width = 50.0
SigT = 1.0
SigS = 0.1 # 0.1, 0.5, 0.99
SigS = 0.5
SigS = 0.99
src  = 1.0



methods = ['step','diamond','stepCharacteristic','linearDiscontinuous']
methods = ['linearDiscontinuous']

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
    
numCells = [90,100,200,400,500,1000,2000,5000]
numCells = [90,100,200,400,700,1000]

sumErrors = getConvergence(numCells,'step',SigT,SigS,src,width)
print(sumErrors)
sumErrors = getConvergence(numCells,'diamond',SigT,SigS,src,width)
print(sumErrors)
sumErrors = getConvergence(numCells,'stepCharacteristic',SigT,SigS,src,width)
print(sumErrors)
sumErrors = getConvergence(numCells,'linearDiscontinuous',SigT,SigS,src,width)
print(sumErrors)
plt.xlabel('# Cells used')
plt.ylabel('% Error Relative to very fine run')
plt.legend(loc='best')

plt.title('S'+str(S.N)+' method for '+str(width)+'cm slab, spatial convergence for scatter = '+str(SigS))

    
plt.show()

"""
plt.title('S'+str(S.N)+' method for '+str(width)+'cm slab, spatial convergence for '+methods[0])
#plt.title('S'+str(S.N)+' method for '+str(width)+'cm slab, using '+str(numCells)+' cells')
plt.xlabel('distance (cm)')
plt.ylabel('avg. scalar flux')
plt.legend(loc='best')
plt.show()



"""









