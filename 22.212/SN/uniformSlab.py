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
	

S = GaussLegendre(2)
#S = GaussLegendre(4)

numCells = 80
width = 20.0
dx = width/numCells
SigT = 1.0
SigS = 0.1 # 0.1, 0.5, 0.99
src  = 1.0

cells = []
for i in range(numCells):
    cells.append(cell(i,SigT,SigS,src,dx))



methods = ['step','diamond','stepCharacteristic','linearDiscontinuous']
for method in methods:
    print(method)
    numIter = 13
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
        plt.plot(x,phi,label=method)

        diff = [oldPhi[i]-phi[i] for i in range(numCells)]
        if max(diff) < 1e-6: 
            converged = True
            print(counter)
        oldPhi = phi[:]

    

plt.title('S'+str(S.N)+' method for '+str(width)+'cm slab, using '+str(numCells)+' cells')
plt.xlabel('distance (cm)')
plt.ylabel('avg. scalar flux')
#plt.legend(loc='best')
plt.show()












