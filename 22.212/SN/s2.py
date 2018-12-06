import matplotlib.pyplot as plt
import numpy as np

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


class slab:
    def __init__(self,cells,width):
        self.cells = cells
        self.width = width
        self.dx = width/len(cells)
    def __str__(slab):
        return "This slab has width "+str(slab.width)+"cm and "+\
                str(len(slab.cells))+" cells with spacing "+str(slab.dx)+"cm"
    def whichCellAmIIn(slab,x):
        for cell in self.cells:
            if cell.L < x < cell.R:
                return cell



##############################################################################
# Define quadrature capabilities
##############################################################################
from numpy.polynomial.legendre import leggauss
class GaussLegendre:
    def __init__(self, N):
        self.N = N
        self.mu, self.wgt = leggauss(N)
	

s2 = GaussLegendre(2)

width = 50.0
numCells = 1000
spacing = width/numCells
SigT = 1.0
SigS = 0.1
src  = 1.0

cellVec = []
for cell_i in range(numCells):
    L = cell_i*spacing
    R = L + spacing
    cellVec.append(cell(cell_i,SigT,SigS,src,spacing))

slab = slab(cellVec,width)




def getQ(cell):
    return cell.src + 0.5*cell.SigS*cell.phi

def diamondDifferenceForward(mu,cell,psiOld,dx):
    denom = (2.0*mu + cell.SigT*dx)
    Q = getQ(cell)
    return ( (2.0*mu - cell.SigT) * psiOld + 2.0*Q*dx ) / denom

def diamondDifferenceBackward(mu,cell,psiOld,dx):
    denom = (2.0*mu - cell.SigT*dx)
    Q = getQ(cell)
    return ( (2.0*mu + cell.SigT) * psiOld - 2.0*Q*dx ) / denom



"""
dx = slab.dx
muF = s2.mu[0]
muB = s2.mu[1]
forwardPsiL  = 0.0
backwardPsiR = 0.0
slabPsi_forward  = [0.0]*(numCells+1)
slabPsi_backward = [0.0]*(numCells+1)

for i in range(numCells):
    cell = slab.cells[i]

    forward_index = i
    backward_index = numCells-i-1

    # calculate psi i+1/2 assuming you know psi i-1/2 for mu > 0
    forwardPsiR  = diamondDifferenceForward( muF,cell,forwardPsiL, dx)

    # calculate psi i-1/2 assuming you know psi i+1/2 for mu < 0
    backwardPsiL = diamondDifferenceBackward(muB,cell,backwardPsiR,dx)


    # This is a bit redundant (lot of overwriting) but it works.
    # slabPsi_forward[0] = psi 1/2  (left)
    # slabPsi_forward[1] = psi 3/2  (right)
    # 
    # slabPsi_forward[1] = psi 3/2  (left)
    # slabPsi_forward[2] = psi 5/2  (right)
    # 
    # ... 
    slabPsi_forward[forward_index] = forwardPsiL
    slabPsi_forward[forward_index+1] = forwardPsiR
    slabPsi_backward[backward_index+1] = backwardPsiR
    slabPsi_backward[backward_index] = backwardPsiL

    # every cell has a psi i-1/2 value and an psi i+1/2 value
    slab.cells[forward_index].PsiL_forward = forwardPsiL
    slab.cells[forward_index].PsiR_forward = forwardPsiR
    # (do this for forward and backward)
    slab.cells[backward_index].PsiL_backward = backwardPsiL
    slab.cells[backward_index].PsiR_backward = backwardPsiR

    forwardPsiL  = forwardPsiR
    backwardPsiR = backwardPsiL 


tau = 0.5


node_x = [slab.cells[0].L]
for cell in slab.cells:
    node_x.append(cell.R)


#plt.plot(node_x,slabPsi_forward)
#plt.plot(node_x,slabPsi_backward)

cell_x = []
psiAvgF = []
psiAvgB = []
for i in range(numCells):
    cell = slab.cells[i]
    cell_x.append(cell.M)
    psiAvgF.append((cell.PsiL_forward+cell.PsiR_forward)*0.5)
    psiAvgB.append((cell.PsiL_backward+cell.PsiR_backward)*0.5)


totalPsi = []
for cell in slab.cells:
    totalPsi.append(psiAvgF[cell.ID]+psiAvgB[cell.ID])

#plt.plot(cell_x,psiAvgF)
#plt.plot(cell_x,psiAvgB)
plt.plot(cell_x,totalPsi)
plt.show()
         




"""

# Forward
psi_in = 0.0
mu = abs(s2.mu[0])
psi_F = []
for i,cell in enumerate(slab.cells[:-1]):
    dx = slab.dx
    Q = getQ(cell) 
    psi_out = ( psi_in * (2*mu - cell.SigT*dx) + 2.0*Q*dx) / \
              ( 2*mu + cell.SigT*dx )
    slab.cells[i+1].psi = psi_out
    psi_F.append(psi_out)
    psi_in = psi_out

# Backward
psi_in = 0.0
mu = abs(s2.mu[1])
psi_B = []
for i,cell in enumerate(slab.cells[:-1]):
    n = len(slab.cells)-i-2
    dx = slab.dx
    Q = getQ(cell) 
    psi_out = ( psi_in * (2*mu - cell.SigT*dx) + 2.0*Q*dx) / \
              ( 2*mu + cell.SigT*dx )
    slab.cells[n].psi = psi_out
    psi_B.append(psi_out)
    psi_in = psi_out


psi_total = psi_F + [psi_B[-i-1] for i in range(len(psi_B))]
plt.plot(psi_total)

plt.show()
















