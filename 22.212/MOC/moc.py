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

def doesItCrossCircle(c,intersections,r):
    A = r.sin*r.sin+r.cos*r.cos
    B = 2.0*(r.x0*r.cos - c.x*r.cos + r.y0*r.sin - c.y*r.sin)
    C = r.x0**2 + c.x*c.x + r.y0**2 + c.y**2 - c.r*c.r

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
        x0 = (cell.R.x-cell.L.x)*(random.random()-0.5)
        y0 = (cell.U.y-cell.D.y)*(random.random()-0.5)
        cos, sin = 2.0*random.random()-1.0, 2.0*random.random()-1.0
    
        r = ray(x0,y0,sin,cos,300.0) # should be 3m length
        r = ray(x0,y0,sin,cos,3.0) # should be 3m length
        
        xPlanes, yPlanes = [cell.L,cell.R], [cell.U,cell.D]
        circles = cell.C

        # figure out where x0 y0 is, then make psi = q[that]/sigmaT
        whereRayIs = whichCircleIsRayIn(r,circles)
        if not whereRayIs: r.psi = sim.q[0]/(4.0*pi*cell.mat.SigmaT)
        else             : r.psi = sim.q[1]/(4.0*pi*cell.C[0].mat.SigmaT)
        rayTracks[rayNum][0].append(["initial whereRayIs",whereRayIs])
        recentDist = [] 
        counter = 0
        while(r.l > 0):
            r.xMax, r.yMax = r.x0+r.cos*r.l, r.y0+r.sin*r.l
            intersections = findAllIntersections(xPlanes,yPlanes,circles,r)
            if intersections == []: print("NO MORE INTERSECTIONS"); return
            firstIntersection = findFirstIntersection(intersections,r)

            whereRayIs = whichCircleIsRayIn(r,circles)
            if not whereRayIs: whereRayIs = cell
             
            if (firstIntersection["t"] > r.l): r.x = r.xMax; r.y = r.yMax; break
            if (firstIntersection["t"] > r.l): r.x = r.xMax; r.y = r.yMax; break
            if weAreStuck(recentDist,firstIntersection): print("GOT STUCK"); return
    
            currentCircle = whichCircleIsRayIn(r,circles)
    
            #if (whereRayIs): plotRaySegment(r,whereRayIs.color,firstIntersection)
            #if (currentCircle): plotRaySegment(r,currentCircle.color,firstIntersection)
            #else: plotRaySegment(r,'#1e4877',firstIntersection)
            updateRay(r,firstIntersection)
               
            length = firstIntersection["dist to get to x1y1"]/np.cos(r.azimuthal)
                
            mat = whereRayIs.mat; idVal = whereRayIs.id
    
            rayTracks[rayNum].append({"length":length,"mat":mat,"idVal":idVal})
    
        
            deltaPsi = (r.psi - (sim.q[idVal]/(4.0*pi*mat.SigmaT)))*(1.0-exp(-mat.SigmaT*length))
            r.psi -= deltaPsi
            #if ( r.l < 250.0 ):
            sim.phi[idVal] += 4.0*pi*deltaPsi
            #print("Delta Psi",deltaPsi)
            #print('%0.03f'%deltaPsi,'%0.03f'%sim.phi[idVal],'%0.03f'%length)
            counter += 1

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
            deltaPsi = (psi - (sim.q[idVal]/(4.0*pi*mat.SigmaT)))*(1.0-exp(-mat.SigmaT*length))
            psi  -= deltaPsi
            #if ( r.l < 250.0 ):
            sim.phi[idVal] += 4.0*pi*deltaPsi
 
    

# Initialize Geometry
radii = [0.5,0.4,0.3,0.2,0.1]
color = ['#ff3366','#cc5980','#998099','#66a6b3','#33cccc']

radii = [0.5,0.3]
color = ['#cc5980','#66a6b3']
radii = [0.39128]
color = ['#cc5980']

# Initialize Materials
#        Sigma T, nuSigma F, Sigma S
mod  = material(3.60, 0.00, 2.10)
fuel = material(0.93, 1.98, 0.30)


rayTracks = []



# box(box_x,box_y,sideLength,radii vec, mod, fuel, color vec for circles, color)
cell = box(0.0,0.0,1.26,radii,mod,fuel,color,'#1e4877')
allRegions = [cell]+cell.C
volumes = [0.0]*len(allRegions)
"""
volumes[0] = cell.sideLength - pi*cell.C[0].r**2
volumes[1] = pi*cell.C[0].r**2 - pi*cell.C[1].r**2
volumes[2] = pi*cell.C[1].r**2 - pi*cell.C[2].r**2
volumes[3] = pi*cell.C[2].r**2 - pi*cell.C[3].r**2
volumes[4] = pi*cell.C[3].r**2 - pi*cell.C[4].r**2
volumes[5] = pi*cell.C[4].r**2
"""
"""
volumes[0] = cell.sideLength - pi*cell.C[0].r**2
volumes[1] = pi*cell.C[0].r**2 - pi*cell.C[1].r**2
volumes[2] = pi*cell.C[1].r**2
"""
volumes[0] = cell.sideLength - pi*cell.C[0].r**2
volumes[1] = pi*cell.C[0].r**2 


k_guess = 1.0

q_guess = [0.0]*len(allRegions)
phi_guess = [0.0]*len(allRegions)
oldFissionSource = [0.0]*len(allRegions)
newFissionSource = [0.0]*len(allRegions)
for i in range(len(q_guess)):
    mat = allRegions[i].mat
    q_guess[i] = ( mat.SigmaS + mat.SigmaF )
    oldFissionSource[i] = mat.SigmaF

converged = False
numRaysPerRun = 1

counter = 0
kVals = []
firstIteration = True
while not converged:
    random.seed(5)
    sim = simulation(phi_guess,q_guess,k_guess) 
    runRays(sim,numRaysPerRun,firstIteration)

    inv_trackLength = 1.0/(300.0*numRaysPerRun)
    sim.phi = [inv_trackLength * x for x in sim.phi]

    phi_Old = sim.phi[:]

    # Update phi
    for i in range(len(sim.q)):
        mat = allRegions[i].mat
        idVal = allRegions[i].id
        sim.phi[i] = sim.phi[i]/(mat.SigmaT*volumes[idVal]) + sim.q[i]/mat.SigmaT

    # Update q
    for i in range(len(q_guess)):
        mat = allRegions[i].mat
        sim.q[i] = ( mat.SigmaS*sim.phi[i] + mat.SigmaF*sim.phi[i] )
        newFissionSource[i] = (1.0/sim.k)*mat.SigmaF*sim.phi[i]

    sim.k = sum(newFissionSource)/sum(oldFissionSource)

    print()
    print(sim.k)
    print(" ".join('%0.5f' % item for item in newFissionSource))
    #print(" ".join('%0.5f' % item for item in sim.phi))
    #print(" ".join('%0.5f' % item for item in q_guess))

    k_guess = sim.k
    kVals.append(sim.k)
    phi_guess = [0.0]*len(allRegions)
    q_guess = sim.q[:]
    # Check if converged
    counter += 1
    if counter > 1: converged = True
    firstIteration = False

#print(rayTracks)
#plotBox(ax,cell)

#plt.xlim(-1,1)
#plt.ylim(-1,1)
#plt.show()
#plt.plot(kVals)
#plt.show()
