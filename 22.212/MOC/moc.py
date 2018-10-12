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
    def __init__(self,phi,k):
        self.phi = phi; self.k = k; self.psi = 1.0 

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
                x_int = r.x0 + t*r.cos
                y_int = r.y0 + t*r.sin
                intersections.append({"xInt":x_int,"yInt":y_int,
                                      "t":t,"bc":None,"dir":None})
    return intersections


def findFirstIntersection(intersections,r):
    t_min = intersections[0]["t"]
    first_intersection = intersections[0]
    for intersection in intersections:
        t = intersection["t"]
        if ( t < t_min ):
            first_intersection = intersection
            t_min = t

    r.x1 = first_intersection["xInt"]
    r.y1 = first_intersection["yInt"]
    dist_traveled = ((r.x1-r.x0)**2+(r.y1-r.y0)**2)**0.5
    first_intersection["dist to get to x1y1"] = dist_traveled
    first_intersection["t"] = t_min
    return first_intersection

def updateRay(r,intersection,counter):
    continue_forward = True
    # do we bounce off ?
    if intersection["bc"] == "ref":
        continue_forward = False
        if   intersection["dir"] == "x": r.cos = -r.cos
        elif intersection["dir"] == "y": r.sin = -r.sin
   
    dist_traveled = intersection["dist to get to x1y1"]

    # Update position and distance to travel
    r.x0 = r.x1
    r.y0 = r.y1
    r.l = r.l - dist_traveled

    return continue_forward


def whichCircleIsRayIn(r,circles):
    potential_circles = []
    for circle in circles:
        if (isRayInCircle(r,circle)):
            potential_circles.append(circle)
    if (len(potential_circles)==0):
        return None
    whereRayIs = potential_circles[0]
    for circle in potential_circles:
        if (circle.r < whereRayIs.r):
            whereRayIs = circle
    return whereRayIs


def isRayInCircle(r,circle):
    smudgeX = r.x0 + 0.001*r.cos
    smudgeY = r.y0 + 0.001*r.sin
    return (smudgeX-circle.x)**2+(smudgeY-circle.y)**2 <= (circle.r)**2


def weAreStuck(recentDist,firstIntersection):
    if len(recentDist) >= 5:
        recentDist.pop(0)
    recentDist.append(firstIntersection["dist to get to x1y1"])


    if sum(recentDist) < 1e-3:
        return True
    return False



def runRays(init_phi,k,iterNum):

    sim = simulation(init_phi,k) 

    # Initialize Ray 
    x0 = (cell.R.x-cell.L.x)*(random.random()-0.5)
    y0 = (cell.U.y-cell.D.y)*(random.random()-0.5)
    cos = 2.0*random.random()-1.0
    sin = 2.0*random.random()-1.0

    r = ray(x0,y0,sin,cos,4.0) # should be 4 length
    
    xPlanes = [cell.L,cell.R]
    yPlanes = [cell.U,cell.D]
    circles = cell.C

    counter = 0
    recentDist = [] 

    while(r.l > 0):
        r.xMax, r.yMax = r.x0+r.cos*r.l, r.y0+r.sin*r.l
    
        intersections = findAllIntersections(xPlanes,yPlanes,circles,r)
        if intersections == []: break

        firstIntersection = findFirstIntersection(intersections,r)
        whereRayIs = whichCircleIsRayIn(r,circles)
        if not whereRayIs: whereRayIs = cell
        if (whereRayIs): plotRaySegment(r,whereRayIs.color,firstIntersection)
        #else:               plotRaySegment(r,'#1e4877',          firstIntersection)

         
        if (firstIntersection["t"] > r.l): r.x = r.xMax; r.y = r.yMax; break
        if weAreStuck(recentDist,firstIntersection): break

        updateRay(r,firstIntersection,counter)
       
        length = firstIntersection["dist to get to x1y1"]/r.azimuthal
        
        mat = whereRayIs.mat; idVal = whereRayIs.id

        incomingPsi = sim.psi

        Q = 0.25/pi * (mat.SigmaS * sim.phi[idVal] + (1.0/sim.k)*mat.SigmaF*sim.phi[idVal]) 
        deltaPsi = (incomingPsi - Q/mat.SigmaT) * (1.0 - exp(-mat.SigmaT*length))
        sim.phi[idVal] += deltaPsi
        if ( r.l < 350.0 ):
            k_numer = 0
            k_denom = 0
            for i in range(len(sim.phi)):
                k_numer += allRegions[i].mat.SigmaF * sim.phi[i] 
                k_denom += allRegions[i].mat.SigmaA * sim.phi[i] 
            sim.k = k_numer * sim.k / k_denom

            sizePhi_inv = 1.0/sum(sim.phi)
            for i in range(len(sim.phi)):
                sim.phi[i] *= sizePhi_inv

            sim.psi = sim.psi - deltaPsi
        
        counter += 1

    return sim.k,sim.phi


    

# Initialize Geometry
radii = [0.1,0.2,0.3,0.4]
color = ['#ee3e32','#f68838','#fbb021','#1b8a5a']

radii = [0.1,0.2,0.3,0.4,0.5]
color = ['#33cccc','#66a6b3','#998099','#cc5980','#ff3366']

# Initialize Materials
#        Sigma T, Sigma A, nuSigma F, Sigma S
fuel = material(50.0, 5.0, 7.0, 28.0)
mod  = material(50.0, 5.0, 7.0, 28.0)

# box(box_x,box_y,sideLength,radii vec, mod, fuel, color vec for circles, color)
cell = box(0.0,0.0,1.0,radii,mod,fuel,color,'#1e4877')
allRegions = [cell]+cell.C



#for i in range(20):
k = 1.0
k_recent = k
numIter = 4 
phi_Init = [1.0]*len(allRegions)
for i in range(numIter):
    #print(i)
    k_recent, phi_recent = runRays(phi_Init,k_recent,i)
    k += k_recent
    phi_Init = phi_recent
k /= numIter

print(k)

plotBox(ax,cell)

plt.xlim(-1,1)
plt.ylim(-1,1)
plt.show()

