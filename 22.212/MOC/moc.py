import matplotlib.pyplot as plt
import random
from geom import *
from plotting import *
from materials import *
from math import pi
from math import exp
import time

random.seed(3)
fig, ax = plt.subplots() 

class simulation():
    def __init__(self,phi,q,k):
        self.phi = phi; self.k = k; self.psi = 100.0; self.q = q

class ray:
    def __init__(self,x,y,sin,cos,length):
        self.x0 = x; self.y0 = y; self.sin = sin; self.cos = cos;
        self.l = length; self.azimuthal = random.random()*pi - pi*0.5;

def findAllIntersections(xPlanes,yPlanes,circles,r):
    intersections = []

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

    for c in circles:
        A = r.sin*r.sin+r.cos*r.cos
        B = 2.0*(r.x0*r.cos - c.x*r.cos + r.y0*r.sin - c.y*r.sin)
        C = r.x0**2 + c.x*c.x + r.y0**2 + c.y**2 - c.r*c.r

        plus_minus = [1.0,-1.0]
        for pm in plus_minus:
            t = (-B + pm*(B*B-4*A*C)**0.5)/(2.0*A)
            if (abs(t.imag) < 1.0e-20 and t.real > 0.0):
                t = t.real
                xInt, yInt = r.x0 + t*r.cos, r.y0 + t*r.sin
                intersections.append({"xInt":xInt, "yInt":yInt,
                                      "t":t, "bc":None, "dir":None})
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
    continue_forward = True
    # do we bounce off ?
    if intersection["bc"] == "ref":
        continue_forward = False
        if   intersection["dir"] == "x": r.cos = -r.cos
        elif intersection["dir"] == "y": r.sin = -r.sin
   
    dist_traveled = intersection["dist to get to x1y1"]

    # Update position and distance to travel
    r.x0, r.y0 = r.x1, r.y1
    r.l = r.l - dist_traveled

    return continue_forward


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
    smudgeX, smudgeY = r.x0 + 0.001*r.cos, r.y0 + 0.001*r.sin
    return (smudgeX-circle.x)**2+(smudgeY-circle.y)**2 <= (circle.r)**2


def weAreStuck(recentDist,firstIntersection):
    if len(recentDist) >= 5: recentDist.pop(0)
    recentDist.append(firstIntersection["dist to get to x1y1"])
    return True if (sum(recentDist) < 1e-3) else False



def runRays(sim,numRaysPerRun):
    for ray in range(numRaysPerRun):
        runRay(sim)


def runRay(sim):
    # Initialize Ray 
    x0 = (cell.R.x-cell.L.x)*(random.random()-0.5)
    y0 = (cell.U.y-cell.D.y)*(random.random()-0.5)
    cos, sin = 2.0*random.random()-1.0, 2.0*random.random()-1.0

    r = ray(x0,y0,sin,cos,4.0) # should be 4 length
    
    xPlanes, yPlanes = [cell.L,cell.R], [cell.U,cell.D]
    circles = cell.C

    recentDist = [] 

    while(r.l > 0):
        r.xMax, r.yMax = r.x0+r.cos*r.l, r.y0+r.sin*r.l
    
        intersections = findAllIntersections(xPlanes,yPlanes,circles,r)
        if intersections == []: break

        firstIntersection = findFirstIntersection(intersections,r)
        whereRayIs = whichCircleIsRayIn(r,circles)
        if not whereRayIs: whereRayIs = cell
        if (whereRayIs): plotRaySegment(r,whereRayIs.color,firstIntersection)
         
        if (firstIntersection["t"] > r.l): r.x = r.xMax; r.y = r.yMax; break
        if weAreStuck(recentDist,firstIntersection): return

        updateRay(r,firstIntersection)
       
        # Account for the fact that we're actually moving in 3d and I'm just 
        # projecting things into 2d
        length = firstIntersection["dist to get to x1y1"]/r.azimuthal
        
        mat = whereRayIs.mat; idVal = whereRayIs.id

        incomingPsi = sim.psi

        deltaPsi = (incomingPsi - sim.q[idVal]/mat.SigmaT)*(1.0 - exp(-mat.SigmaT*length))
        sim.psi = sim.psi - deltaPsi
        if ( r.l < 350.0 ):
            sim.phi[idVal] += deltaPsi/(mat.SigmaT*volumes[idVal])
        

    

# Initialize Geometry
radii = [0.1,0.2,0.3,0.4]
color = ['#ee3e32','#f68838','#fbb021','#1b8a5a']

radii = [0.1,0.2,0.3,0.4,0.5]
color = ['#33cccc','#66a6b3','#998099','#cc5980','#ff3366']

# Initialize Materials
#        Sigma T, Sigma A, nuSigma F, Sigma S
fuel = material(50.0, 7.0, 2.4*5.0, 28.0)
mod  = material(50.0, 7.0, 2.4*5.0, 28.0)

fuel = material(0.93, 0.00, 0.98, 1.00)
mod  = material(3.60, 00.0, 0.00, 3.30)



# box(box_x,box_y,sideLength,radii vec, mod, fuel, color vec for circles, color)
cell = box(0.0,0.0,1.0,radii,mod,fuel,color,'#1e4877')
allRegions = [cell]+cell.C
volumes = [0.0]*len(allRegions)
volumes[0] = cell.sideLength - pi*cell.C[4].r**2
volumes[1] = pi*cell.C[4].r**2 - pi*cell.C[3].r**2
volumes[2] = pi*cell.C[3].r**2 - pi*cell.C[2].r**2
volumes[3] = pi*cell.C[2].r**2 - pi*cell.C[1].r**2
volumes[4] = pi*cell.C[1].r**2 - pi*cell.C[0].r**2
volumes[5] = pi*cell.C[0].r**2


k_guess = 1.0
q_guess = [2.0/(4.0*pi)]*1 + [0.35/(4.0*pi)]*5
phi_guess = [q_guess[0]/mod.SigmaT]*1 + [q_guess[1]/fuel.SigmaT]*5

converged = False
numRaysPerRun = 20

counter = 0
while not converged:
    sim = simulation(phi_guess,q_guess,k_guess) 
    runRays(sim,numRaysPerRun)
    
    # Update q
    sim.q[0] = 1.0/(4.0*pi) * (  (1.0/sim.k) * allRegions[0].mat.SigmaF * sim.phi[0] + 
                             (allRegions[0].mat.SigmaS * sim.phi[0] ) )
    sim.q[1] = 1.0/(4.0*pi) * (  (1.0/sim.k) * allRegions[1].mat.SigmaF * sim.phi[1] + 
                             (allRegions[1].mat.SigmaS * sim.phi[1] ) )
    sim.q[2] = 1.0/(4.0*pi) * (  (1.0/sim.k) * allRegions[2].mat.SigmaF * sim.phi[2] + 
                             (allRegions[2].mat.SigmaS * sim.phi[2] ) )
    sim.q[3] = 1.0/(4.0*pi) * (  (1.0/sim.k) * allRegions[3].mat.SigmaF * sim.phi[3] + 
                             (allRegions[3].mat.SigmaS * sim.phi[3] ) )
    sim.q[4] = 1.0/(4.0*pi) * (  (1.0/sim.k) * allRegions[4].mat.SigmaF * sim.phi[4] + 
                             (allRegions[4].mat.SigmaS * sim.phi[4] ) )
    sim.q[5] = 1.0/(4.0*pi) * (  (1.0/sim.k) * allRegions[5].mat.SigmaF * sim.phi[5] + 
                             (allRegions[5].mat.SigmaS * sim.phi[5] ) )

    sim.phi[0] = sim.q[0]/allRegions[0].mat.SigmaT
    sim.phi[1] = sim.q[1]/allRegions[1].mat.SigmaT
    sim.phi[2] = sim.q[2]/allRegions[2].mat.SigmaT
    sim.phi[3] = sim.q[3]/allRegions[3].mat.SigmaT
    sim.phi[4] = sim.q[4]/allRegions[4].mat.SigmaT
    sim.phi[5] = sim.q[5]/allRegions[5].mat.SigmaT

    # Normalize Phi
    inv_sum_phi = 1.0 / sum(sim.phi)
    sim.phi = [inv_sum_phi * x for x in sim.phi]


    # Check if converged
    # converged = True
    counter += 1
    if counter > 10: converged = True

print(sim.phi)

plotBox(ax,cell)
plt.xlim(-1,1)
plt.ylim(-1,1)
