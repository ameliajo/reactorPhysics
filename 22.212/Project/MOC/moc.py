import matplotlib.pyplot as plt
import random
from geom import *
from plotting import *
from materials import *
from math import pi
from math import exp
import numpy as np
from colors import *
import time
import numpy as np
import sys


numRaysPerRun = 100
rayDist = 200.0
deadZone = 20.0
fissionSourceError = 0.01
kError = 0.1
print("Running",numRaysPerRun,"rays for a distance of",rayDist,"with a deadzone of",deadZone)

fig, ax = plt.subplots() 

class ray:
    def __init__(self,x,y,length):
        self.x0 = x; self.y0 = y; self.l = length; self.psi = None
        polar = 2.0*pi*random.random()
        self.cos_azi = random.random()
        self.sin_azi = np.sin(np.arccos(self.cos_azi))
        self.sin, self.cos = self.sin_azi*np.sin(polar), self.sin_azi*np.cos(polar)


def initToZero(nGroups,nCellRegs,nCells):
    return [[[0.0 for k in range(nGroups)] for j in range(nCellRegs)] for i in range(nCells)]

def getCell(cells,r):
    """
    Since we're working with a 3x3 lattice, we need to find which of the 9 we're in. 
    Input:  Vector filled with cell objects, and our ray object
    Output: Correct cell object, List of X-Planes, List of Y-Planes, and list 
            of circles in the cell 
    """
    cell, counter = None, 0
    while not cell:
        for cell in cells:
            if r.x0 >= cell.L.x and r.x0 <= cell.R.x and r.y0 >= cell.D.y and r.y0 <= cell.U.y:
                return cell, [cell.L,cell.R], [cell.U,cell.D], cell.C

        # smudge factor bc corners were out to get me
        smudge = 1.0e-5 * (-1)**counter * (1+counter)
        r.x0, r.y0 = r.x0+smudge*r.cos, r.y0+smudge*r.sin
        counter += 1



def doesItCrossCircle(c,intersections,r):
    """
    Update the intersections vector to account for whether the ray passes 
    through given circle c. I want to see if my ray is going to cross the circle 
    in question. So I take (x-xc)^2 + (y-yc)^2 = r^2, parameterize 
    x = x0 + t*cos, y = y0 + t*sin, and solving for t will require solving the 
    quadratic formula below. If t is  real and positive, then it's in our path 
    and we should add it to the intersections list.

    Inputs: circle object, existing list of intersections, ray object
    Output: amended list of intersections
    """
    A = r.sin**2 + r.cos**2
    B = 2*(r.x0*r.cos - c.x*r.cos + r.y0*r.sin - c.y*r.sin)
    C = r.x0**2 - 2*r.x0*c.x + c.x*c.x + r.y0**2 - 2*r.y0*c.y + c.y**2 - c.r**2

    for pm in [1.0,-1.0]:
        if B*B - 4*A*C < 0: continue
        t = (-B + pm*(B*B-4*A*C)**0.5)/(2*A)
        if (abs(t.imag) < 1.0e-20 and t.real > 0.0):
            t, xInt, yInt = t.real, r.x0 + t*r.cos, r.y0 + t*r.sin
            intersections.append({"x":xInt,"y":yInt,"t":t,"bc":None,"dir":None})
    return intersections



def findAllIntersections(xPlanes,yPlanes,circles,r):
    """
    Find all potential intersections that ray r will go through. Consider three
    cases: ray is outside a circle, but inside cell, ray is in a circle but not
    the innermost circle, and ray is in the innermost circle. 1st case means 
    that we are only interested in wall intersections and intersection with the
    largest circle. 2nd case means we're interested in the boundary of our 
    current circle, and the boundary of the next smaller ring. 3rd case means
    we only need to check the boundary of our current circle.
    """
    intersections = []
    whereRayIs = whichCircleIsRayIn(r,circles)
    if not whereRayIs:
        for xPlane in xPlanes:
            t = (xPlane.x - r.x0)/r.cos
            if t <= 0: continue
            xInt, yInt = round(r.x0+t*r.cos,10), round(r.y0+t*r.sin,10)
            intersections.append({"x":xInt,"y":yInt,"t":t,"bc":xPlane.bc,"dir":"x"})

        for yPlane in yPlanes:
            t = (yPlane.y - r.y0)/r.sin
            if t <= 0: continue
            xInt, yInt = round(r.x0+t*r.cos,10), round(r.y0+t*r.sin,10)
            intersections.append({"x":xInt,"y":yInt,"t":t,"bc":yPlane.bc,"dir":"y"})

        return doesItCrossCircle(circles[0],intersections,r)
    
    # If I'm on the innermost circle, I can only cross this inner circle to go out
    if whereRayIs.id == len(circles):
        return doesItCrossCircle(whereRayIs,intersections,r)
    
    # If I'm on any other circle, I can either cross my boundary or the one inside me
    for c in [circles[whereRayIs.id-1],circles[whereRayIs.id]]:
        intersections = doesItCrossCircle(c,intersections,r)

    return intersections


def findFirstIntersection(xPlanes,yPlanes,circles,r,recentDist):
    """
    Find the first intersection that ray r is going to encounter. Uses the 
    function ``findAllIntersections'' to get all possible/logical candidates.
    We look at the intersection with the minimum t (parametrization variable),
    which should indicate which it hits first.
    """
    intersections = findAllIntersections(xPlanes,yPlanes,circles,r)
    if intersections == []: print("NO MORE INTERSECTIONS"); return None
    times = [i["t"] for i in intersections]
    firstInt = intersections[times.index(min(times))]
    r.x1, r.y1 = firstInt["x"], firstInt["y"]
    firstInt["dist"] = ((r.x1-r.x0)**2+(r.y1-r.y0)**2)**0.5/r.sin_azi
    if weAreStuck(recentDist,firstInt): print("GOT STUCK"); return None
    return firstInt,recentDist
 


def updateRay(r,intersection):
    """
    Update ray by moving it to the next boundary, and potentially changing its
    direction (if we hit a reflective boundary).
    """
    if intersection["bc"] == "ref":     # Do we bounce off?
        if   intersection["dir"] == "x": r.cos = -r.cos 
        elif intersection["dir"] == "y": r.sin = -r.sin
   
    r.x0, r.y0 = r.x1 + 1e-7*r.cos, r.y1 + 1e-7*r.sin
    r.l -= intersection["dist"]



def whichCircleIsRayIn(r,circles):
    outsideRings = [c for c in circles if (r.x0-c.x)**2+(r.y0-c.y)**2 <= c.r**2]
    if outsideRings == []: return None
    radiusList   = [c.r for c in outsideRings]
    return outsideRings[radiusList.index(min(radiusList))]


def weAreStuck(recentDist,firstIntersect):
    if len(recentDist) < 5 : return False
    recentDist.pop(0); recentDist.append(firstIntersect["dist"])
    return True if (sum(recentDist) < 1e-8) else False


def runRays(numRaysPerRun,k,phi,Q):
    distInFuel, distInMod = 0.0, 0.0
    for ray in range(numRaysPerRun):
        dFuel,dMod = runRay(ray,k,phi,Q)
        distInFuel += dFuel; 
        distInMod  += dMod; 

    return distInMod,distInFuel     


def runRay(rayNum,k,phi,Q):

    # Initialize Ray 
        
    x0 = (cells[2].R.x-cells[0].L.x)*(random.random()) + cells[0].L.x
    y0 = (cells[8].U.y-cells[0].D.y)*(random.random()) + cells[0].D.y
    
    r = ray(x0,y0,rayDist) # should be 3m length

    cell,xPlanes,yPlanes,circles = getCell(cells,r)

    recentDist = [] 
    distInFuel,distInMod = 0.0, 0.0

    while(r.l > 0):

        cell,xPlanes,yPlanes,circles = getCell(cells,r)
        r.xMax, r.yMax = r.x0+r.cos*r.l, r.y0+r.sin*r.l

        firstIntersect,recentDist = findFirstIntersection(xPlanes,yPlanes,circles,r,recentDist)
        if not firstIntersect: print('something bad'); return distInFuel,distInMod

        whereRayIs = whichCircleIsRayIn(r,circles)
        if whereRayIs: idVal = whereRayIs.id
        else:          idVal = 0; whereRayIs = cell; 
        mat = whereRayIs.mat

        if not r.psi: 
            r.psi = [Q[cell.id][idVal][g] for g in range(nGroups)]
             
        if (firstIntersect["t"] > r.l): r.x = r.xMax; r.y = r.yMax; break
    
        if len(sys.argv) > 1:
            if sys.argv[1] == "diagram" and whereRayIs:
                plotRaySegment(r,whereRayIs.color,firstIntersect)

        updateRay(r,firstIntersect)

        deltaPsiVec = []
        for g in range(nGroups):
            deltaPsi_g = (r.psi[g] - Q[cell.id][idVal][g]) * \
                         (1.0 - exp(-mat.SigT[g]*firstIntersect["dist"]))



            r.psi[g] = r.psi[g] - deltaPsi_g
            deltaPsiVec.append(deltaPsi_g)

        if ( r.l < (rayDist-deadZone) ):
            for g in range(nGroups):
                phi[cell.id][idVal][g] += 4.0*pi*deltaPsiVec[g]
                        
            if idVal != 0: distInFuel += firstIntersect["dist"]
            else:          distInMod  += firstIntersect["dist"]


    return distInFuel,distInMod



def runMOC(oldFissSrc,newFissSrc,cells,k,phi,Q):

    converged = False; 
    kVals = []; counter = 0; invDist = 0.0
    invDist = 0.0

    while not converged:
        random.seed(1)
        distInMod,distInFuel = runRays(numRaysPerRun,k,phi,Q)
        invDist = 1.0/(distInMod+distInFuel)


        phi = [[[invDist*phi[i][j][g] for g in range(nGroups)]  \
                   for j in range(nCellRegs)] for i in range(nCells) ]

        
        # Update phi and normalize phi
        kNumer, kDenom = 0.0, 0.0
        for i in range(nCells):
            for j in range(nCellRegs):
                for g in range(nGroups):
                    mat = allReg[i][j].mat
                    V = (distInMod)*invDist if j == 0 else (distInFuel)*invDist
                    phi[i][j][g] = phi[i][j][g]/(mat.SigT[g]*V) + Q[i][j][g]*4.0*pi

                    newFissSrc[i][j][g] = mat.SigF[g]*phi[i][j][g]
                    #newFissSrc[i][j][g] = getFissionIntoG(g,mat.chi,mat.SigF,phi[i][j])/k

                    kNumer += newFissSrc[i][j][g]
                    kDenom += oldFissSrc[i][j][g]
 

        diff_k = abs(k - kNumer/kDenom)
        k = kNumer / kDenom
 
        phiSum = sum([sum([sum(i) for i in j]) for j in phi])
        phi = [[[groupTerm/phiSum for groupTerm in regTerm] for regTerm in cellTerm] for cellTerm in phi]

        totalDiff = 0.0
        for i in range(nCells):
            for j in range(nCellRegs):
                for g in range(nGroups):
                    mat = allReg[i][j].mat

                    newFissSrc[i][j][g] = mat.SigF[g]*phi[i][j][g]/k
                    #newFissSrc[i][j][g] = getFissionIntoG(g,mat.chi,mat.SigF,phi[i][j])/k


                    sigS = getScatteringIntoG(g,mat.SigS_matrix,phi[i][j])
                    Q[i][j][g] = (sigS+newFissSrc[i][j][g]) / (mat.SigT[g]*4*pi)

                    totalDiff += abs(oldFissSrc[i][j][g] - newFissSrc[i][j][g])
                    oldFissSrc[i][j][g] = newFissSrc[i][j][g]

    
        print("------  Run # ",counter,"       k-eff",k,"        Error in Fission Source",totalDiff)
        kVals.append(k)
    
        # Check if converged
        if totalDiff < fissionSourceError and diff_k < kError: 
            for i in phi: print("MOD  ",[float("%0.2E"%k) for k in i[0]])
            for i in phi: print("FUEL ",[float("%0.2E"%k) for k in i[1]])

            f = open("fluxMOC.py","w+")
            for i in range(9):
                f.write("MOC_modFlux"+str(i)+"  = "+str([float("%.8f"%flux) for flux in phi[i][0]])+"\n")
                f.write("MOC_fuelFlux"+str(i)+" = "+str([float("%.8f"%flux) for flux in phi[i][1]])+"\n")
                f.write("\n\n")
            f.close()


            return kVals,phi

        phi = initToZero(nGroups,nCellRegs,nCells)
        counter += 1

    






##############################################################################
# Initialize Geometry and Materials
##############################################################################
#radii = [0.39128,0.3]
radii = [0.39128]
sideLen = 1.26

xPlanes = [xPlane(0.0,'ref'),xPlane(sideLen,'vac'),xPlane(2.0*sideLen,'vac'),xPlane(3.0*sideLen,'ref')]
yPlanes = [yPlane(0.0,'ref'),yPlane(sideLen,'vac'),yPlane(2.0*sideLen,'vac'),yPlane(3.0*sideLen,'ref')]


# Store all pincell objects in vector
i = 0
cells = []
for y in range(3):
    for x in range(3):
        cells.append(box(sideLen*x+sideLen*0.5+xPlanes[0].x,               \
                         sideLen*y+sideLen*0.5+yPlanes[0].y,radii,modV[i], \
                         fuelV[i],i,x,y,xPlanes,yPlanes))
        i += 1

nCellRegs = 1 + len(radii)

allReg = []
for cell in cells:
    allReg.append([cell]+cell.C)


nGroups = len(modV[0].SigT)
nCells  = 9

oldFissSrc = initToZero(nGroups,nCellRegs,nCells)
newFissSrc = initToZero(nGroups,nCellRegs,nCells)

Q          = initToZero(nGroups,nCellRegs,nCells)
phi        = initToZero(nGroups,nCellRegs,nCells)

for c in range(nCells):
    for i in range(nCellRegs):
        mat = allReg[c][i].mat
        for g in range(nGroups):
            sigS = getScatteringIntoG(g,mat.SigS_matrix,[1.0]*nGroups)
            sigF = getFissionIntoG(g,mat.chi,mat.SigF,[1.0]*nGroups)
            sigF = mat.SigF[g]
            oldFissSrc[c][i][g] = sigF
            Q[c][i][g] = ( sigS + sigF )/(mat.SigT[g]*4.0*pi)

Qsum = sum([sum([sum(i) for i in j]) for j in Q])
Q = [[[groupTerm/Qsum for groupTerm in regTerm] for regTerm in cellTerm] for cellTerm in Q]

k = 1.0


##############################################################################
# Actually running MOC
##############################################################################

kVals, phi = runMOC(oldFissSrc,newFissSrc,cells,k,phi,Q)




##############################################################################
# Plotting and Alerting of Negative Flux
##############################################################################

if len(sys.argv) > 1:
    if sys.argv[1] == "keff":
        f1 = plt.figure(1)
        plt.plot(kVals,'ro')
        f1.show()
        input("\nPress Enter to close plot")
    if sys.argv[1] == "diagram":
        for cell in cells:
            plotBox(ax,cell)
        plt.show()

print()
for i in phi:
    for j in i:
        for k in j:
            if k <= 0.0:
                print("Got a negative flux value! :( ",k)

print()
print(kVals[-1])




"""

"""


