import matplotlib.pyplot as plt
import random
from geom import *
from plotting import *
from materials import *
from math import pi
from math import exp
import time
import numpy as np

fig, ax = plt.subplots() 

class simulation():
    def __init__(self,phi,q,k):
        self.phi = phi; self.k = k; self.q = q

class ray:
    def __init__(self,x,y,sin,cos,length):
        self.x0 = x; self.y0 = y; self.sin = sin; self.cos = cos;
        self.l = length; 
        self.azimuthal = random.random()*pi - pi*0.5;
        #self.azimuthal = 0.05


def whichCellAmIIn(r,cells):
    for cell in cells:
        if r.x0 > cell.L.x and r.x0 < cell.R.x and r.y0 > cell.D.y and r.y0 < cell.U.y:
            return cell
    return None

def doesItCrossCircle(c,intersections,r):

    #print(circles[0].x,circles[0].y,circles[0].r,circles[0].id)
    A = r.sin*r.sin+r.cos*r.cos
    B = 2.0*(r.x0*r.cos - c.x*r.cos + r.y0*r.sin - c.y*r.sin)
    C = r.x0**2 - 2.0*r.x0*c.x + c.x*c.x + r.y0**2 -2.0*r.y0*c.y + c.y**2 - c.r*c.r

    for pm in [1.0,-1.0]:
        t = (-B + pm*(B*B-4*A*C)**0.5)/(2.0*A)
        if (abs(t.imag) < 1.0e-20 and t.real > 0.0):
            t = t.real
            xInt, yInt = r.x0 + t*r.cos, r.y0 + t*r.sin
            intersections.append({"xInt":xInt, "yInt":yInt,"t":t, "bc":None, "dir":None})
    return intersections



def findAllIntersections(xPlanes,yPlanes,circles,r):
    intersections = []
    whereRayIs = whichCircleIsRayIn(r,circles)
    if not whereRayIs:
        for xPlane in xPlanes:
            t = (xPlane.x - r.x0)/r.cos
            if ( t > 0 ):
                xInt, yInt = round(r.x0+t*r.cos,10), round(r.y0+t*r.sin,10)
                intersections.append({"xInt":xInt,  "yInt":yInt, "t":t,
                                      "bc":xPlane.bc,"dir":"x"})

        for yPlane in yPlanes:
            t = (yPlane.y - r.y0)/r.sin
            if ( t > 0 ):
                xInt, yInt = round(r.x0+t*r.cos,10), round(r.y0+t*r.sin,10)
                intersections.append({"xInt":xInt,  "yInt":yInt, "t":t,
                                      "bc":yPlane.bc,"dir":"y"})

        return doesItCrossCircle(circles[0],intersections,r)

    
    if whereRayIs.id == len(circles):
        return doesItCrossCircle(whereRayIs,intersections,r)

    
    for c in [circles[whereRayIs.id-1],circles[whereRayIs.id]]:
        intersections = doesItCrossCircle(c,intersections,r)
    return intersections


def findFirstIntersection(intersections,r):
    tMin = intersections[0]["t"]
    firstInt = intersections[0]
    for intersection in intersections:
        t = intersection["t"]
        if ( t < tMin ):
            firstInt = intersection
            tMin = t

    r.x1, r.y1 = firstInt["xInt"], firstInt["yInt"]
    dist_traveled = ((r.x1-r.x0)**2+(r.y1-r.y0)**2)**0.5
    firstInt["t"], firstInt["dist to get to x1y1"] = tMin, dist_traveled
    return firstInt

def updateRay(r,intersection):
    # do we bounce off ?
    if intersection["bc"] == "ref":
        if   intersection["dir"] == "x": r.cos = -r.cos
        elif intersection["dir"] == "y": r.sin = -r.sin
   
    dist_traveled = intersection["dist to get to x1y1"]

    # Update position and distance to travel
    r.x0, r.y0 = r.x1+0.0001*r.cos, r.y1+0.0001*r.sin
    r.l = r.l - dist_traveled



def whichCircleIsRayIn(r,circles):
    potential_circles = []
    for circle in circles:
        if (isRayInCircle(r,circle)): potential_circles.append(circle)
    if ( potential_circles == [] ): return None

    whereRayIs = potential_circles[0]
    for circle in potential_circles:
        if (circle.r < whereRayIs.r): whereRayIs = circle
    return whereRayIs


def isRayInCircle(r,circle):
    return (r.x0-circle.x)**2+(r.y0-circle.y)**2 <= (circle.r)**2


def weAreStuck(recentDist,firstIntersection):
    if len(recentDist) >= 5: recentDist.pop(0)
    recentDist.append(firstIntersection["dist to get to x1y1"])
    return True if (sum(recentDist) < 1e-3) else False



def runRays(sim,numRaysPerRun,firstIteration):
    for ray in range(numRaysPerRun):
        if firstIteration:
            rayTracks.append([["RayNumber"+str(ray)]])
        runRay(sim,ray,firstIteration)


def runRay(sim,rayNum,firstIteration):

    if firstIteration:
        # Initialize Ray 
        x0 = (cell1.R.x-cell1.L.x)*(random.random()-0.5)
        y0 = (cell1.U.y-cell1.D.y)*(random.random()-0.5)
        cos, sin = 2.0*random.random()-1.0, 2.0*random.random()-1.0
    
        r = ray(x0,y0,sin,cos,300.0) # should be 3m length

        cell = whichCellAmIIn(r,cells)
        xPlanes, yPlanes = [cell.L,cell.R], [cell.U,cell.D]
        circles = cell.C

        # figure out where x0 y0 is, then make psi = q[that]/sigmaT
        whereRayIs = whichCircleIsRayIn(r,circles)
        if not whereRayIs: r.psi = sim.q[cell.id][0]/(4.0*pi*cell.mat.SigmaT)
        else             : r.psi = sim.q[cell.id][1]/(4.0*pi*cell.C[0].mat.SigmaT)
        rayTracks[rayNum][0].append(["initial whereRayIs",whereRayIs])
        recentDist = [] 
        while(r.l > 0):
            #print("\nIn Cell",whichCellAmIIn(r,cells).id)
            cell = whichCellAmIIn(r,cells)
            xPlanes, yPlanes = [cell.L,cell.R], [cell.U,cell.D]
            circles = cell.C
            

            r.xMax, r.yMax = r.x0+r.cos*r.l, r.y0+r.sin*r.l
            intersections = findAllIntersections(xPlanes,yPlanes,circles,r)
            if intersections == []: print("NO MORE INTERSECTIONS"); return
            firstIntersection = findFirstIntersection(intersections,r)

            whereRayIs = whichCircleIsRayIn(r,circles)
            if whereRayIs: mat = whereRayIs.mat; idVal = whereRayIs.id
            if not whereRayIs: idVal = 0; whereRayIs = cell
            #print("LocationID",idVal)
             
            if (firstIntersection["t"] > r.l): r.x = r.xMax; r.y = r.yMax; break
            if (firstIntersection["t"] > r.l): r.x = r.xMax; r.y = r.yMax; break
            if weAreStuck(recentDist,firstIntersection): print("GOT STUCK"); return
    
            currentCircle = whichCircleIsRayIn(r,circles)
            if (whereRayIs): plotRaySegment(r,whereRayIs.color,firstIntersection)
            if (currentCircle): plotRaySegment(r,currentCircle.color,firstIntersection)
            else: plotRaySegment(r,'#1e4877',firstIntersection)
            updateRay(r,firstIntersection)
            length = firstIntersection["dist to get to x1y1"]/np.cos(r.azimuthal)
            rayTracks[rayNum].append({"length":length,"mat":mat,"idVal":idVal,"l":r.l})
        
            deltaPsi = (r.psi - (sim.q[cell.id][idVal]/(4.0*pi*mat.SigmaT)))*(1.0-exp(-mat.SigmaT*length))
            r.psi -= deltaPsi
            if ( r.l < 250.0 ):
                sim.phi[cell.id][idVal] += 4.0*pi*deltaPsi

    else:
        tracksForThisRun = rayTracks[rayNum]
        whereRayIs = tracksForThisRun[0][1]
        psi = 0.0;
        if not whereRayIs: psi = sim.q[0]/(4.0*pi*cell.mat.SigmaT)
        else             : psi = sim.q[1]/(4.0*pi*cell.C[0].mat.SigmaT)
        for i in range(len(tracksForThisRun)-1):
            mat    = tracksForThisRun[i+1]["mat"]
            idVal  = tracksForThisRun[i+1]["idVal"]
            length = tracksForThisRun[i+1]["length"]
            l      = tracksForThisRun[i+1]["l"]
            deltaPsi = (psi - (sim.q[idVal]/(4.0*pi*mat.SigmaT)))*(1.0-exp(-mat.SigmaT*length))
            psi  -= deltaPsi
            if ( l < 250.0 ):
                sim.phi[idVal] += 4.0*pi*deltaPsi
 
    

# Initialize Geometry

#radii = [0.5,0.4,0.3,0.2,0.1]
#color = ['#ff3366','#cc5980','#998099','#66a6b3','#33cccc']

radii = [0.5,0.3]
color = ['#cc5980','#66a6b3']

#radii = [0.39128]
#color = ['#cc5980']

# Initialize Materials
#        Sigma T, nuSigma F, Sigma S
mod  = material(3.60, 0.00, 2.10)
fuel = material(0.93, 1.98, 0.30)

xPlanes = [xPlane(-0.63,'ref'),xPlane(0.63,'vac'),xPlane(1.89,'vac'),xPlane(3.15,'ref')]
yPlanes = [yPlane(-0.63,'ref'),yPlane(0.63,'vac'),yPlane(1.89,'vac'),yPlane(3.15,'ref')]

rayTracks = []

# box(box_x,box_y,sideLength,radii vec, mod, fuel, color vec for circles, color)
cell0 = box(0.00,0.00,1.26,radii,mod,fuel,color,'#1e4877',0,xPlanes[0],xPlanes[1],yPlanes[0],yPlanes[1])
cell1 = box(1.26,0.00,1.26,radii,mod,fuel,color,'#1e4877',1,xPlanes[1],xPlanes[2],yPlanes[0],yPlanes[1])
cell2 = box(2.52,0.00,1.26,radii,mod,fuel,color,'#1e4877',2,xPlanes[2],xPlanes[3],yPlanes[0],yPlanes[1])
cell3 = box(0.00,1.26,1.26,radii,mod,fuel,color,'#1e4877',3,xPlanes[0],xPlanes[1],yPlanes[1],yPlanes[2])
cell4 = box(1.26,1.26,1.26,radii,mod,fuel,color,'#1e4877',4,xPlanes[1],xPlanes[2],yPlanes[1],yPlanes[2])
cell5 = box(2.52,1.26,1.26,radii,mod,fuel,color,'#1e4877',5,xPlanes[2],xPlanes[3],yPlanes[1],yPlanes[2])
cell6 = box(0.00,2.52,1.26,radii,mod,fuel,color,'#1e4877',6,xPlanes[0],xPlanes[1],yPlanes[2],yPlanes[3])
cell7 = box(1.26,2.52,1.26,radii,mod,fuel,color,'#1e4877',7,xPlanes[1],xPlanes[2],yPlanes[2],yPlanes[3])
cell8 = box(2.52,2.52,1.26,radii,mod,fuel,color,'#1e4877',8,xPlanes[2],xPlanes[3],yPlanes[2],yPlanes[3])
cells = [cell0,cell1,cell2,cell3,cell4,cell5,cell6,cell7,cell8]
allRegionsCell1 = [cell1]+cell1.C
volumes = [0.0]*len(allRegionsCell1)
"""
volumes[0] = cell1.sideLength - pi*cell1.C[0].r**2
volumes[1] = pi*cell1.C[0].r**2 - pi*cell1.C[1].r**2
volumes[2] = pi*cell1.C[1].r**2 - pi*cell1.C[2].r**2
volumes[3] = pi*cell1.C[2].r**2 - pi*cell1.C[3].r**2
volumes[4] = pi*cell1.C[3].r**2 - pi*cell1.C[4].r**2
volumes[5] = pi*cell1.C[4].r**2
"""
volumes[0] = cell1.sideLength - pi*cell1.C[0].r**2
volumes[1] = pi*cell1.C[0].r**2 - pi*cell1.C[1].r**2
volumes[2] = pi*cell1.C[1].r**2
"""
volumes[0] = cell1.sideLength - pi*cell1.C[0].r**2
volumes[1] = pi*cell1.C[0].r**2 
"""


k_guess = 1.0

q_guess = []
phi_guess = []
oldFissionSource = []
newFissionSource = []
for i in range(len(cells)):
    q_guess.append([0.0]*len(allRegionsCell1))
    phi_guess.append([0.0]*len(allRegionsCell1))
    oldFissionSource.append([0.0]*len(allRegionsCell1))
    newFissionSource.append([0.0]*len(allRegionsCell1))

#q_guess = [0.0]*len(allRegionsCell1)
#phi_guess, oldFissionSource, newFissionSource = q_guess[:], q_guess[:], q_guess[:]
for c in range(len(cells)):
    for i in range(len(allRegionsCell1)):
        mat = allRegionsCell1[i].mat
        q_guess[c][i] = ( mat.SigmaS + mat.SigmaF )
        oldFissionSource[c][i] = mat.SigmaF

converged = False
firstIteration = True
numRaysPerRun = 1 

counter = 0
kVals = []
while not converged:
    random.seed(5)
    sim = simulation(phi_guess,q_guess,k_guess) 
    runRays(sim,numRaysPerRun,firstIteration)

    inv_trackLength = 1.0/(300.0*numRaysPerRun)
    for cell in cells:
        sim.phi[cell.id] = [inv_trackLength * x for x in sim.phi[cell.id]]


    # Update phi
    for cell in cells:
        for i in range(len(allRegionsCell1)):
            mat = allRegionsCell1[i].mat
            idVal = allRegionsCell1[i].id
            sim.phi[cell.id][i] = sim.phi[cell.id][i]/(mat.SigmaT*volumes[idVal]) + sim.q[cell.id][i]/mat.SigmaT
            sim.q[cell.id][i] = ( mat.SigmaS*sim.phi[cell.id][i] + mat.SigmaF*sim.phi[cell.id][i] )
            newFissionSource[cell.id][i] = (1.0/sim.k)*mat.SigmaF*sim.phi[cell.id][i]

    kNumer = 0
    kDenom = 0
    for i in range(len(newFissionSource)):
        kNumer += sum(newFissionSource[i])
        kDenom += sum(oldFissionSource[i])
    sim.k = kNumer / kDenom
    #sim.k = sum(sum(newFissionSource))/sum(sum(oldFissionSource))

    print()
    print(sim.k)
    #print(" ".join('%0.5f' % item for item in newFissionSource))

    k_guess = sim.k
    kVals.append(sim.k)
    phi_guess = [0.0]*len(allRegionsCell1)
    q_guess = sim.q[:]
    # Check if converged
    counter += 1
    if counter > 0: converged = True
    firstIteration = False

#print(rayTracks)
for cell in cells:
    plotBox(ax,cell)

plt.xlim(-1,3.5)
plt.ylim(-1,3.5)
#plt.show()
#plt.plot(kVals)
plt.show()
