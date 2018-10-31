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

fig, ax = plt.subplots() 

class simulation():
    def __init__(self,k,q):
        self.k = k; 

        self.phi = []
        for cell_i in q:
            cellEntry = []
            for location_j in range(len(cell_i)):
                locationEntry = []
                for g in range(nGroups):
                    locationEntry.append(0.0)
                cellEntry.append(locationEntry)
            self.phi.append(cellEntry)


        self.q = []
        for cell_i in q:
            cellEntry = []
            for location_j in range(len(cell_i)):
                locationEntry = []
                for g in range(nGroups):
                    locationEntry.append(q[i][j][g])
                cellEntry.append(locationEntry)
            self.q.append(cellEntry)



class ray:
    def __init__(self,x,y,polar,cos_azi,length):
        self.x0 = x; self.y0 = y; self.l = length; self.psi = None
        self.cos_azi = cos_azi
        sin_azi = np.sin(np.arccos(cos_azi))
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
        #cell = whichCellAmIIn(r,cells)
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
    firstInt["dist"] = ((r.x1-r.x0)**2+(r.y1-r.y0)**2)**0.5/r.cos_azi
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
    radiusList = [c.r for c in outsideRings]
    return outsideRings[radiusList.index(min(radiusList))]


def weAreStuck(recentDist,firstIntersect):
    if len(recentDist) < 5 : return False
    recentDist.pop(0)
    recentDist.append(firstIntersect["dist"])
    return True if (sum(recentDist) < 1e-8) else False


def runRays(sim,numRaysPerRun,firstIteration=False):
    distInFuel, distInMod = 0.0, 0.0
    for ray in range(numRaysPerRun):
        if firstIteration: rayTracks.append([["RayNumber"+str(ray)]])
        dFuel,dMod = runRay(sim,ray,firstIteration)
        if firstIteration:
            distInFuel += dFuel; distInMod  += dMod; 

    return distInMod,distInFuel     


def runRay(sim,rayNum,firstIteration):

    if firstIteration:
        # Initialize Ray 
        
        x0 = (cells[2].R.x-cells[0].L.x)*(random.random()) + cells[0].L.x
        y0 = (cells[8].U.y-cells[0].D.y)*(random.random()) + cells[0].D.y
        polar = 2.0*pi*random.random()
        cos_azimuthal = random.random()
        
        r = ray(x0,y0,polar,cos_azimuthal,30.0) # should be 3m length

        cell,xPlanes,yPlanes,circles = getCell(cells,r)


        recentDist = [] 
        distInFuel,distInMod = 0.0, 0.0
        psi = []
        while(r.l > 0):

            cell,xPlanes,yPlanes,circles = getCell(cells,r)
            
            r.xMax, r.yMax = r.x0+r.cos*r.l, r.y0+r.sin*r.l

            firstIntersect,recentDist = findFirstIntersection(xPlanes,yPlanes,circles,r,recentDist)
            if not firstIntersect: print('something bad'); return distInFuel,distInMod

            whereRayIs = whichCircleIsRayIn(r,circles)
            if whereRayIs: idVal = whereRayIs.id
            else:          idVal = 0; whereRayIs = cell; 
            mat = whereRayIs.mat
            if not psi: 
                rayTracks[rayNum][0].append(["initial whereRayIs",whereRayIs])
                psi = []
                for g in range(nGroups):
                    psi.append(sim.q[cell.id][idVal][g]/(mat.SigT[g]))
             
            if (firstIntersect["t"] > r.l): r.x = r.xMax; r.y = r.yMax; break
    
            #if (whereRayIs): plotRaySegment(r,whereRayIs.color,firstIntersect)

            updateRay(r,firstIntersect)

            rayTracks[rayNum].append({"length":firstIntersect["dist"],"mat":mat,
                                      "idVal":idVal,"l":r.l,"cellID":cell.id})
        
            deltaPsiVec = []
            for g in range(nGroups):
                deltaPsi_g = (psi[g] - (sim.q[cell.id][idVal][g]/(mat.SigT[g]))) *  \
                             (1.0 - exp(-mat.SigT[g]*firstIntersect["dist"]))
                psi[g] = psi[g] - deltaPsi_g
                deltaPsiVec.append(deltaPsi_g)

            if ( r.l < 250.0 ):
                for g in range(nGroups):
                    sim.phi[cell.id][idVal][g] += 4.0*pi*deltaPsiVec[g]
                if whereRayIs.id != 0: distInFuel += firstIntersect["dist"]
                else:                  distInMod  += firstIntersect["dist"]


        rayTracks[rayNum][0] += [["distInFuel",distInFuel],["distInMod",distInMod]]

        return distInFuel,distInMod

    else:
        tracks = rayTracks[rayNum]
        mat, cellId = rayTracks[0][2]["mat"], rayTracks[0][2]["cellID"]
        whereRayIs = tracks[0][1]

        psi = []
        for g in range(nGroups):
            psi.append( sim.q[cellId][1][g]/(mat.SigT[g]) if whereRayIs else   \
                        sim.q[cellId][0][g]/(mat.SigT[g]) )



        for i in range(len(tracks)-1):
            mat    = tracks[i+1]["mat"   ]; idVal  = tracks[i+1]["idVal"]
            length = tracks[i+1]["length"]; l      = tracks[i+1]["l"    ]
            cellID = tracks[i+1]["cellID"]; 
            deltaPsiVec = []
            for g in range(nGroups):
                deltaPsiVec.append( (psi[g] - (sim.q[cellID][idVal][g]/(mat.SigT[g])))*    \
                           (1.0 - exp(-mat.SigT[g]*length)) )
                psi[g] -= deltaPsiVec[g]

            if ( l < 250.0 ):
                for g in range(nGroups):
                    sim.phi[cellID][idVal][g] += 4.0*pi*deltaPsiVec[g]

        distInFuel, distInMod = tracks[0][2][1], tracks[0][3][1]
        return distInFuel,distInMod
 
    

# Initialize Geometry

#radii = [0.5,0.4,0.3,0.2,0.1]
#radii = [0.4,0.3,0.2,0.1]
radii = [0.39128,0.3]
#radii = [0.39128,0.3]
#radii = [0.39128]

cellColors, circleColors = findColors(radii)


# Initialize Materials
#        Sigma T, nuSigma F, Sigma S

#mod  = material([3.60], [0.00], [2.10])
#fuel = material([0.93], [1.98], [0.30])

#mod  = material([2.50], [0.00], [2.10])
#fuel = material([1.12], [0.98], [0.19])

#mod  = material([2.50,2.50,2.50], [0.00,0.00,0.00], [2.10,2.10,2.10])
#fuel = material([1.12,1.12,1.12], [0.98,0.98,0.98], [0.19,0.19,0.19])


S_matrix_M = [ [ 0.6, 0.5 ], [ 0.4, 0.3 ] ]
S_matrix_F = [ [ 0.1, 0.2 ], [ 0.3, 0.4 ] ]

S_matrix_M = [ [ 0.21, 0.02, 0.23 ], [ 0.06, 0.95, 0.94 ], [ 0.13, 0.12, 0.11 ] ]
S_matrix_F = [ [ 0.19, 0.18, 0.17 ], [ 0.16, 0.15, 0.14 ], [ 0.13, 0.12, 0.11 ] ]

chi_M = [ 0.00, 0.00, 0.00 ]
chi_F = [ 0.05, 0.05, 0.90 ]

mod  = material([2.50,2.50,2.50], [0.00,0.00,0.00], [2.10,2.10,2.10],S_matrix_M,chi_M)
fuel = material([1.12,1.12,1.12], [0.98,0.98,0.98], [0.19,0.19,0.19],S_matrix_F,chi_F)

nGroups = len(mod.SigT)



xPlanes = [xPlane(-0.63,'ref'),xPlane(0.63,'vac'),xPlane(1.89,'vac'),xPlane(3.15,'ref')]
yPlanes = [yPlane(-0.63,'ref'),yPlane(0.63,'vac'),yPlane(1.89,'vac'),yPlane(3.15,'ref')]

rayTracks = []

# box(box_x,box_y,sideLength,radii vec, mod, fuel, color vec for circles, color)
sideLength = 1.26
counter = 0
cells = []
for y in range(3):
    for x in range(3):
        cell = box(1.26*x,1.26*y,radii,mod,fuel,circleColors[counter],cellColors[counter],counter,xPlanes[x],xPlanes[x+1],yPlanes[y],yPlanes[y+1])
        counter += 1
        cells.append(cell)

allRegInCell = [cells[0]]+cells[0].C
numRegionsPerCell = 1 + len(radii)

volumes = getVolumes(numRegionsPerCell,cells[0].C,sideLength)

k_guess = 1.0

q_guess = []
oldFissionSource = []
newFissionSource = []

for i in range(len(cells)):
    cellEntry = []
    for j in range(numRegionsPerCell):
        regionEntry = []
        for k in range(nGroups):
            regionEntry.append(0.0)
        cellEntry.append(regionEntry)
    q_guess.append(cellEntry[:])

for i in range(len(cells)):
    cellEntry = []
    for j in range(numRegionsPerCell):
        regionEntry = []
        for k in range(nGroups):
            regionEntry.append(0.0)
        cellEntry.append(regionEntry)
    oldFissionSource.append(cellEntry[:])

for i in range(len(cells)):
    cellEntry = []
    for j in range(numRegionsPerCell):
        regionEntry = []
        for k in range(nGroups):
            regionEntry.append(0.0)
        cellEntry.append(regionEntry)
    newFissionSource.append(cellEntry[:])

for c in range(len(cells)):
    for i in range(numRegionsPerCell):
        mat = allRegInCell[i].mat
        for g in range(nGroups):
            fakeFlux = [1.0]*nGroups
            sigS = getScatteringIntoG(g,mat.SigS_matrix,fakeFlux)
            sigF = getFissionIntoG(g,mat.chi,mat.SigF,fakeFlux)
            q_guess[c][i][g] = ( sigS + sigF )/(4.0*pi)
            oldFissionSource[c][i][g] = sigF



converged = False
firstIteration = True
numRaysPerRun = 100

random.seed(1)
counter = 0
inv_trackLength = 0.0
kVals = []
print()
#print("Q GUESS")
#print(q_guess)
sim = simulation(k_guess,q_guess) 
#print()
#print("Q")
#print(sim.q)
#print()


while not converged:

    if counter == 0: 
        distInMod,distInFuel = runRays(sim,numRaysPerRun,True)
        inv_trackLength = 1.0/(distInMod+distInFuel)
    else:                      runRays(sim,numRaysPerRun)

    for i,sim_phi_cell in enumerate(sim.phi):
        for j,x in enumerate(sim_phi_cell):
            for g in range(nGroups):
                sim.phi[i][j][g] = inv_trackLength*x[g]

    # Update phi
    for cell in cells:
        for i in range(numRegionsPerCell):
            for g in range(nGroups):
                mat, V = allRegInCell[i].mat, volumes[allRegInCell[i].id]
                sim.phi[cell.id][i][g] = sim.phi[cell.id][i][g]/(mat.SigT[g]*V) + sim.q[cell.id][i][g]*4.0*pi/mat.SigT[g]
                sigS = getScatteringIntoG(g,mat.SigS_matrix,sim.phi[cell.id][i])
                sim.q[cell.id][i][g] = (sigS + mat.SigF[g]*sim.phi[cell.id][i][g])/(4.0*pi)
                #sim.q[cell.id][i][g] = (mat.SigS[g]*sim.phi[cell.id][i][g] + mat.SigF[g]*sim.phi[cell.id][i][g])/(4.0*pi)
                newFissionSource[cell.id][i][g] = mat.SigF[g]*sim.phi[cell.id][i][g]


    kNumer, kDenom = 0.0, 0.0
    for i in range(len(newFissionSource)):
        for j in range(len(newFissionSource[0])):
            for g in range(nGroups):
                kNumer += newFissionSource[i][j][g]
                kDenom += oldFissionSource[i][j][g]
    sim.k = kNumer / kDenom

    for cell in cells:
        for i in range(numRegionsPerCell):
            for g in range(nGroups):
                oldFissionSource[cell.id][i][g] = newFissionSource[cell.id][i][g]/sim.k
                sim.q[cell.id][i][g] /= sim.k 

    #print("Iteration #: ",counter,"   k",sim.k)

    kVals.append(sim.k)

    sim.phi = []
    for cell_i in sim.q:
        cellEntry = []
        for location_j in range(len(cell_i)):
            locationEntry = []
            for g in range(nGroups):
                locationEntry.append(0.0)
            cellEntry.append(locationEntry)
        sim.phi.append(cellEntry)




    # Check if converged
    counter += 1
    if counter > 20: converged = True
    firstIteration = False

#correctK = [1.0393502572032347, 1.0356661397728641, 1.0325926732412258, 1.0301681484434255, 1.0282467337804195, 1.026718950054905, 1.0255013959039376, 1.0245297679016676, 1.0237539963860283, 1.0231347788746885, 1.0226410604093656]
#wrong = False
#for i in range(len(kVals)):
#    if abs(kVals[i]-correctK[i]) > 1e-10:
#        print("\nOH NO NOT GOOD\n")
#        wrong = True
#if not wrong: print("\nWhew! That was a close one\n")
#




#print(kVals)
#print(rayTracks)
#for cell in cells:
#    plotBox(ax,cell)
#plt.show()
#plt.plot(kVals)
#plt.show()

f1 = plt.figure(1)
plt.plot(kVals,'ro')
f1.show()
input()
"""
"""




