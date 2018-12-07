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
        self.psiL  = 0.0
        self.psiR  = 0.0
        self.width = width
        self.L     = ID*width
        self.M     = (ID+0.5)*width
        self.R     = (ID+1)*width
    def __str__(cell):
        return "Cell #"+str(cell.ID)+" and spans "+str(cell.L)+" - "+str(cell.R)




##############################################################################
# Define quadrature capabilities
##############################################################################
from numpy.polynomial.legendre import leggauss
class GaussLegendre:
    def __init__(self, N):
        self.N = N
        self.mu, self.wgt = leggauss(N)
	

s2 = GaussLegendre(2)

numCells = 300
width = 50.0
dx = width/numCells
SigT = 1.0
SigS = 0.1 # 0.1, 0.5, 0.99
src  = 1.0

cells = []
for i in range(numCells):
    L = i * dx
    R = L + dx 
    cells.append(cell(i,SigT,SigS,src,dx))



methods = ['step','diamond','step-characteristic']
for method in methods:
    numIter = 1
    for i in range(numIter):
        broom(s2,cells,dx,method)
        phi = []
        x = []
        for cell in cells:
            phi.append(cell.phi)
            x.append(cell.M)
        invMaxVal = 1.0/max(phi)
        phi = [x*invMaxVal for x in phi]
        plt.plot(x,phi,label=method)
    
plt.title('Sn methods for '+str(width)+'cm slab, using '+str(numCells)+' cells')
plt.xlabel('distance (cm)')
plt.ylabel('avg. scalar flux')
plt.legend(loc='best')
plt.show()












