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


numRaysPerRun = 250
rayDist = 100.0
deadZone = 10.0
fissionSourceError = 0.01
print("Running",numRaysPerRun,"rays for a distance of",rayDist,"with a deadzone of",deadZone)

fig, ax = plt.subplots() 

class simulation():
    def __init__(self,k,q):
        self.k = k; 
        self.phi = [[[       0.0 for g in range(nGroups)] for j in range(nCellRegs)] \
                                                          for i in range(nCells) ]
        self.q   = [[[q[i][j][g] for g in range(nGroups)] for j in range(nCellRegs)] \
                                                          for i in range(nCells) ]


class ray:
    def __init__(self,x,y,length):
        self.x0 = x; self.y0 = y; self.l = length; self.psi = None
        polar = 2.0*pi*random.random()
        self.cos_azi = random.random()
        sin_azi = np.sin(np.arccos(self.cos_azi))
        self.sin_azi = sin_azi
        self.sin, self.cos = sin_azi*np.sin(polar), sin_azi*np.cos(polar)


def getCell(cells,r):
    """
    Since we're working with a 3x3 lattice, we need to find which of these 9 we
    are in. 
    Input:  Vector filled with cell objects, and our ray object
    Output: Correct cell object, List of X-Planes, List of Y-Planes, and list 
            of circles in the cell 
    """
    cell, counter = None, 0
    while not cell:
        for cell in cells:
            if r.x0 >= cell.L.x and r.x0 <= cell.R.x and r.y0 >= cell.D.y and r.y0 <= cell.U.y:
                return cell, [cell.L,cell.R], [cell.U,cell.D], cell.C

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
    #firstInt["dist"] = ((r.x1-r.x0)**2+(r.y1-r.y0)**2)**0.5/r.cos_azi
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
    outsideRings = [c   for c in circles if (r.x0-c.x)**2+(r.y0-c.y)**2 <= c.r**2]
    if outsideRings == []: return None
    radiusList   = [c.r for c in outsideRings]
    return outsideRings[radiusList.index(min(radiusList))]


def weAreStuck(recentDist,firstIntersect):
    if len(recentDist) < 5 : return False
    recentDist.pop(0); recentDist.append(firstIntersect["dist"])
    return True if (sum(recentDist) < 1e-8) else False


def runRays(sim,numRaysPerRun,rayTracks,firstIteration=False):
    distInFuel, distInMod = 0.0, 0.0
    for ray in range(numRaysPerRun):
        if firstIteration: rayTracks.append([["RayNumber"+str(ray)]])
        dFuel,dMod = runRay(sim,ray,firstIteration,rayTracks)
        if firstIteration:
            distInFuel += dFuel; distInMod  += dMod; 

    return distInMod,distInFuel     


def runRay(sim,rayNum,firstIteration,rayTracks):

    if firstIteration:
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
                rayTracks[rayNum][0].append(["initial whereRayIs",whereRayIs])
                r.psi = []
                for g in range(nGroups):
                    r.psi.append(sim.q[cell.id][idVal][g]/(mat.SigT[g]))
             
            if (firstIntersect["t"] > r.l): r.x = r.xMax; r.y = r.yMax; break
    
            if len(sys.argv) > 1:
                if sys.argv[1] == "diagram" and whereRayIs:
                    plotRaySegment(r,whereRayIs.color,firstIntersect)

            updateRay(r,firstIntersect)

            rayTracks[rayNum].append({"length":firstIntersect["dist"],"mat":mat,
                                      "idVal":idVal,"l":r.l,"cellID":cell.id})
        
            deltaPsiVec = []
            for g in range(nGroups):
                deltaPsi_g = (r.psi[g] - (sim.q[cell.id][idVal][g]/(mat.SigT[g]))) *  \
                             (1.0 - exp(-mat.SigT[g]*firstIntersect["dist"]))
                r.psi[g] = r.psi[g] - deltaPsi_g
                deltaPsiVec.append(deltaPsi_g)

            if ( r.l < (rayDist-deadZone) ):
                for g in range(nGroups):
                    sim.phi[cell.id][idVal][g] += 4.0*pi*deltaPsiVec[g]
                inFuel = False
                for sigFVal in whereRayIs.mat.SigF:
                    if abs(sigFVal) > 1e-12:
                        inFuel = True
                        
                #if whereRayIs.id != 0: print("This should be fuel",whereRayIs.mat.SigF)
                #if inFuel: print("This should be fuel",whereRayIs.mat.SigF)
                if inFuel: distInFuel += firstIntersect["dist"]
                else:      distInMod  += firstIntersect["dist"]


                #if whereRayIs.id != 0: distInFuel += firstIntersect["dist"]
                #else:                  distInMod  += firstIntersect["dist"]


        rayTracks[rayNum][0] += [["distInFuel",distInFuel],["distInMod",distInMod]]

        return distInFuel,distInMod

    else:
        tracks = rayTracks[rayNum]
        mat, cellId = rayTracks[0][2]["mat"], rayTracks[0][2]["cellID"]
        whereRayIs = tracks[0][1]

        psi = [ sim.q[cellId][1][g]/(mat.SigT[g]) if whereRayIs else      \
                sim.q[cellId][0][g]/(mat.SigT[g]) for g in range(nGroups) ]


        for i in range(len(tracks)-1):
            mat    = tracks[i+1]["mat"   ]; idVal  = tracks[i+1]["idVal"]
            length = tracks[i+1]["length"]; l      = tracks[i+1]["l"    ]
            ID     = tracks[i+1]["cellID"]; 
            deltaPsiVec = []
            
            for g in range(nGroups):
                deltaPsiVec.append( (psi[g] - (sim.q[ID][idVal][g]/(mat.SigT[g])))*    \
                           (1.0 - exp(-mat.SigT[g]*length)) )
                psi[g] -= deltaPsiVec[g]

            if ( l < (rayDist-deadZone) ):
                for g in range(nGroups):
                    sim.phi[ID][idVal][g] += 4.0*pi*deltaPsiVec[g]

        distInFuel, distInMod = tracks[0][2][1], tracks[0][3][1]
        return distInFuel,distInMod
 

def runMOC(sim,k_guess,q_guess,oldFissSrc,newFissSrc,cells,volumes):

    converged = False; firstIteration = True
    rayTracks = []; kVals = []; counter = 0; invDist = 0.0
    invDist = 0.0
    goodDistInMod = 0.0
    goodDistInFuel = 0.0

    while not converged:

        random.seed(1)
        distInMod,distInFuel = runRays(sim,numRaysPerRun,rayTracks,firstIteration)
        if firstIteration: invDist = 1.0/(distInMod+distInFuel)
        if firstIteration: goodDistInMod = distInMod
        if firstIteration: goodDistInFuel = distInFuel

        sim.phi = [[[invDist*sim.phi[i][j][g] for g in range(nGroups)]  \
                   for j in range(nCellRegs)] for i in range(nCells) ]

        
        # Update phi and normalize phi
        for cell in cells:
            for j in range(nCellRegs):
                for g in range(nGroups):
                    mat, V = allRegInCell[j].mat, volumes[allRegInCell[j].id]
                    V2 = (goodDistInMod)*invDist if j == 0 else (goodDistInFuel)*invDist
                    V3 = (goodDistInMod)*invDist if j != 0 else (goodDistInFuel)*invDist
                    #print(V/sum(volumes),V2,V3,j) 
                    sim.phi[cell.id][j][g] = sim.phi[cell.id][j][g]/(mat.SigT[g]*V) \
                                           + sim.q[cell.id][j][g]*4.0*pi/mat.SigT[g]

        kNumer, kDenom = 0.0, 0.0
        for i in range(nCells):
            for j in range(nCellRegs):
                for g in range(nGroups):
                    mat = allRegInCell[j].mat
                    newFissSrc[cell.id][j][g] = mat.SigF[g]*sim.phi[cell.id][j][g]
                    kNumer += newFissSrc[i][j][g]
                    kDenom += oldFissSrc[i][j][g]
        
        sim.k = kNumer / kDenom
 
        phiSum = sum([sum([sum(i) for i in j]) for j in sim.phi])
        for cell in cells:
            for j in range(nCellRegs):
                for g in range(nGroups):
                    sim.phi[cell.id][j][g] /= phiSum
                    sigS = getScatteringIntoG(g,mat.SigS_matrix,sim.phi[cell.id][j])
                    sim.q[cell.id][j][g] = (sigS+(1.0/sim.k)*mat.SigF[g]*sim.phi[cell.id][j][g]) / (4*pi)
                    newFissSrc[cell.id][j][g] = (1.0/sim.k)*mat.SigF[g]*sim.phi[cell.id][j][g]

        """
        kNumer, kDenom = 0.0, 0.0
        for i in range(nCells):
            for j in range(nCellRegs):
                for g in range(nGroups):
                    mat, V = allRegInCell[j].mat, volumes[allRegInCell[j].id]
                    V2 = (goodDistInMod)*invDist if j == 0 else (goodDistInFuel)*invDist
                    kNumer += (V*mat.SigF[g]*sim.phi[i][j][g])
                    kDenom += (V*mat.SigA[g]*sim.phi[i][j][g])
    
        sim.k = kNumer / kDenom
        """
        
        totalDiff = 0.0
        for c in cells:
            for i in range(nCellRegs):
                for g in range(nGroups):
                    totalDiff += abs(oldFissSrc[c.id][i][g] - newFissSrc[c.id][i][g]/sim.k)
                    oldFissSrc[c.id][i][g] = newFissSrc[c.id][i][g]/sim.k
                    sim.q[cell.id][j][g] = (sigS + mat.SigF[g]/sim.k*sim.phi[cell.id][j][g])/(4*pi)

    
        print("------  Run # ",counter,"       k-eff",sim.k,"        Error in Fission Source",totalDiff)
        kVals.append(sim.k)
    
        # Check if converged
        if totalDiff < fissionSourceError: 
            converged = True
            #print(sim.phi)
            for i in sim.phi:
                print("MOD  ",[float("%0.2E"%k) for k in i[0]])
            print()
            for i in sim.phi:
                print("FUEL ",[float("%0.2E"%k) for k in i[1]])
        else: sim.phi = [[[0.0 for g in range(nGroups)] for j in range(nCellRegs)] \
                                                        for i in range(nCells)]

        counter += 1
        # This tells me to not actually trace the rays again, but rather rely on my
        # existing / stored information
        firstIteration = False
        

    return kVals
    

##############################################################################
# Initialize Geometry and Materials
##############################################################################
radii = [0.39128]
# Side length of cell
sideLen = 1.26

xPlanes = [xPlane(-0.63,'ref'),xPlane(0.63,'vac'),xPlane(1.89,'vac'),xPlane(3.15,'ref')]
yPlanes = [yPlane(-0.63,'ref'),yPlane(0.63,'vac'),yPlane(1.89,'vac'),yPlane(3.15,'ref')]

# Store all pincell objects in vector
i = 0
cells = []
for y in range(3):
    for x in range(3):
        cells.append(box(sideLen*x,sideLen*y,radii,modV[i],fuelV[i],i,x,y,xPlanes,yPlanes))
        i += 1

# Since I'm assuming that my cells are all the same, I just use the material
# info for my first cell to represent the rest. 
allRegInCell = [cells[0]]+cells[0].C
nCellRegs = 1 + len(radii)

volumes = getVolumes(nCellRegs,cells[0].C,sideLen)

nGroups = len(modV[0].SigT)
nCells  = 9

q_guess    = [[[0.0 for k in range(nGroups)] for j in range(nCellRegs)]  \
                    for i in range(nCells)]
oldFissSrc = [[[0.0 for k in range(nGroups)] for j in range(nCellRegs)]  \
                    for i in range(nCells)]
newFissSrc = [[[0.0 for k in range(nGroups)] for j in range(nCellRegs)]  \
                    for i in range(nCells)]

for c in range(nCells):
    for i in range(nCellRegs):
        mat = allRegInCell[i].mat
        for g in range(nGroups):
            sigS = getScatteringIntoG(g,mat.SigS_matrix,[1.0]*nGroups)
            sigF = getFissionIntoG(g,mat.chi,mat.SigF,[1.0]*nGroups)
            q_guess[c][i][g] = ( sigS + sigF )/(4.0*pi)
            oldFissSrc[c][i][g] = sigF



##############################################################################
# Actually running MOC
##############################################################################

k_guess = 1.0
sim = simulation(k_guess,q_guess) 

kVals = runMOC(sim,k_guess,q_guess,oldFissSrc,newFissSrc,cells,volumes)




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
for i in sim.phi:
    for j in i:
        for k in j:
            if k <= 0.0:
                print("Got a negative flux value! :( ",k)

print()
print(kVals[-1])



