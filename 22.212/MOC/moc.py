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
    def __init__(self,phi,q,k):
        self.phi = phi; self.k = k; self.q = q

class ray:
    def __init__(self,x,y,sin,cos,azimuthal,length):
        self.x0 = x; self.y0 = y; self.sin = sin; self.cos = cos;
        self.l = length; 
        #self.azimuthal = random.random()*pi - pi*0.5;
        self.azimuthal = azimuthal

def getCell(cells,r):
    cell = None
    counter = 0
    #print(r.x0,r.y0)
    while not cell:
        cell = whichCellAmIIn(r,cells)
        if not cell:
            r.x0 += (1.0e-5*r.cos)*(-1)**counter*(1+counter)
            r.y0 += (1.0e-5*r.sin)*(-1)**counter*(1+counter)
            counter += 1
            #plt.show()
            #break

    xPlanes, yPlanes = [cell.L,cell.R], [cell.U,cell.D]
    circles = cell.C
    return cell,xPlanes,yPlanes,circles


def whichCellAmIIn(r,cells):
    for cell in cells:
        if r.x0 >= cell.L.x and r.x0 <= cell.R.x and r.y0 >= cell.D.y and r.y0 <= cell.U.y:
            return cell
    return None


def doesItCrossCircle(c,intersections,r):
    A = r.sin*r.sin+r.cos*r.cos
    B = 2.0*(r.x0*r.cos - c.x*r.cos + r.y0*r.sin - c.y*r.sin)
    C = r.x0**2 - 2.0*r.x0*c.x + c.x*c.x + r.y0**2 -2.0*r.y0*c.y + c.y**2 - c.r*c.r
    #print(A,B,C)

    for pm in [1.0,-1.0]:
        if (B*B-4*A*C < 0): continue
        t = (-B + pm*(B*B-4*A*C)**0.5)/(2.0*A)
        if (abs(t.imag) < 1.0e-20 and t.real > 0.0):
            t = t.real
            xInt, yInt = r.x0 + t*r.cos, r.y0 + t*r.sin
            intersections.append({"xInt":xInt,"yInt":yInt,"t":t,
                                  "bc":None,"dir":None})
    return intersections



def findAllIntersections(xPlanes,yPlanes,circles,r):
    intersections = []
    whereRayIs = whichCircleIsRayIn(r,circles)
    if not whereRayIs:
        for xPlane in xPlanes:
            t = (xPlane.x - r.x0)/r.cos
            if t <= 0: continue
            xInt, yInt = round(r.x0+t*r.cos,10), round(r.y0+t*r.sin,10)
            intersections.append({"xInt":xInt,"yInt":yInt,"t":t,"bc":xPlane.bc,"dir":"x"})

        for yPlane in yPlanes:
            t = (yPlane.y - r.y0)/r.sin
            if t <= 0: continue
            xInt, yInt = round(r.x0+t*r.cos,10), round(r.y0+t*r.sin,10)
            intersections.append({"xInt":xInt,"yInt":yInt,"t":t,"bc":yPlane.bc,"dir":"y"})

        return doesItCrossCircle(circles[0],intersections,r)
    
    # If I'm on the innermost circle, I can only cross this inner circle to go out
    if whereRayIs.id == len(circles):
        return doesItCrossCircle(whereRayIs,intersections,r)
    
    # If I'm on any other circle, I can either cross my boundary or the one inside me
    for c in [circles[whereRayIs.id-1],circles[whereRayIs.id]]:
        intersections = doesItCrossCircle(c,intersections,r)

    return intersections


def findFirstIntersection(xPlanes,yPlanes,circles,r,recentDist):
    intersections = findAllIntersections(xPlanes,yPlanes,circles,r)
    if intersections == []: print("NO MORE INTERSECTIONS"); return None
    tMin = intersections[0]["t"]
    firstInt = intersections[0]
    for intersection in intersections:
        t = intersection["t"]
        if ( t < tMin ):
            firstInt = intersection
            tMin = t

    r.x1, r.y1 = firstInt["xInt"], firstInt["yInt"]
    dist_traveled = ((r.x1-r.x0)**2+(r.y1-r.y0)**2)**0.5/np.cos(r.azimuthal)
    firstInt["t"], firstInt["dist to get to x1y1"] = tMin, dist_traveled

    if weAreStuck(recentDist,firstInt): print("GOT STUCK"); return None
    return firstInt,recentDist
 


def updateRay(r,intersection):
    if intersection["bc"] == "ref":     # Do we bounce off?
        if   intersection["dir"] == "x": r.cos = -r.cos
        elif intersection["dir"] == "y": r.sin = -r.sin
   
    r.x0, r.y0 = r.x1+1.0e-7*r.cos, r.y1+1.0e-7*r.sin
    r.l = r.l - intersection["dist to get to x1y1"]



def whichCircleIsRayIn(r,circles):
    potential_circles = []
    for circle in circles:
        if (isRayInCircle(r,circle)): potential_circles.append(circle)
    if potential_circles == []: return None

    whereRayIs = potential_circles[0]
    for circle in potential_circles:
        if (circle.r < whereRayIs.r): whereRayIs = circle
    return whereRayIs


def isRayInCircle(r,circle):
    return (r.x0-circle.x)**2+(r.y0-circle.y)**2 <= (circle.r)**2


def weAreStuck(recentDist,firstIntersection):
    if len(recentDist) < 5 :
        return False
    if len(recentDist) >= 5: recentDist.pop(0)
    recentDist.append(firstIntersection["dist to get to x1y1"])
    if (sum(recentDist) < 1e-8): 
        print("We are stuck")
    return True if (sum(recentDist) < 1e-8) else False



def runRays(sim,numRaysPerRun,firstIteration):
    distInFuel = 0.0
    distInMod  = 0.0
    for ray in range(numRaysPerRun):
        #print("Ray Number",ray)
        if firstIteration:
            rayTracks.append([["RayNumber"+str(ray)]])
        dFuel,dMod = runRay(sim,ray,firstIteration)
        distInFuel += dFuel
        distInMod  += dMod
    return distInMod,distInFuel     
    #print("Dist in fuel",distInFuel)
    #print("Dist in mod ",distInMod)


def runRay(sim,rayNum,firstIteration):

    if firstIteration:
        # Initialize Ray 
        x0 = (cells[2].R.x-cells[0].L.x)*(random.random()) + cells[0].L.x
        y0 = (cells[8].U.y-cells[0].D.y)*(random.random()) + cells[0].D.y
        polar = 2.0*pi*random.random()
        azimuthal = pi*random.random()-(pi*0.5)

        #cos, sin = 2.0*random.random()-1.0, 2.0*random.random()-1.0

        #x0 = (cells[0].R.x-cells[0].L.x)*(random.random()) + cells[0].L.x
        #y0 = (cells[0].U.y-cells[0].D.y)*(random.random()) + cells[0].D.y

        r = ray(x0,y0,np.cos(azimuthal)*np.sin(polar),np.cos(azimuthal)*np.cos(polar),azimuthal,300.0) # should be 3m length

        cell,xPlanes,yPlanes,circles = getCell(cells,r)
        #print("just got cells")


        # figure out where x0 y0 is, then make psi = q[that]/sigmaT
        whereRayIs = whichCircleIsRayIn(r,circles)
        if not whereRayIs: r.psi = sim.q[cell.id][0]/(4.0*pi*cell.mat.SigmaT)
        else             : r.psi = sim.q[cell.id][1]/(4.0*pi*cell.C[0].mat.SigmaT)

        #print("found where ray is, going into loop")
        rayTracks[rayNum][0].append(["initial whereRayIs",whereRayIs])
        recentDist = [] 
        distInFuel = 0.0
        distInMod = 0.0
        while(r.l > 0):

            #print("at beginning of loop")
            cell,xPlanes,yPlanes,circles = getCell(cells,r)
            #print("found cells")
            
            r.xMax, r.yMax = r.x0+r.cos*r.l, r.y0+r.sin*r.l

            firstIntersection,recentDist = findFirstIntersection(xPlanes,yPlanes,circles,r,recentDist)
            #print("found first intersection")
            if not firstIntersection: print('something bad'); return

            whereRayIs = whichCircleIsRayIn(r,circles)
            #print("found where ray is")
            if whereRayIs: idVal = whereRayIs.id
            if not whereRayIs: idVal = 0; whereRayIs = cell; 
            mat = whereRayIs.mat
             
            if (firstIntersection["t"] > r.l): r.x = r.xMax; r.y = r.yMax; break
    
            #currentCircle = whichCircleIsRayIn(r,circles)
            #if (whereRayIs): plotRaySegment(r,whereRayIs.color,firstIntersection)
            #if (currentCircle): plotRaySegment(r,currentCircle.color,firstIntersection)
            #else: plotRaySegment(r,cell.color,firstIntersection)

            updateRay(r,firstIntersection)

            #print("updated ray")
            length = firstIntersection["dist to get to x1y1"]#/np.cos(r.azimuthal)

            rayTracks[rayNum].append({"length":length,"mat":mat,"idVal":idVal,"l":r.l,"cellID":cell.id})
            #print("appended to ray tracks")
        
            deltaPsi = (r.psi - (sim.q[cell.id][idVal]/(4.0*pi*mat.SigmaT)))*(1.0-exp(-mat.SigmaT*length))
            #print("calculated delta psi")
            r.psi -= deltaPsi
            if ( r.l < 250.0 ):
                sim.phi[cell.id][idVal] += 4.0*pi*deltaPsi
                if whereRayIs.id == 1:
                    distInFuel += firstIntersection["dist to get to x1y1"]
                else:
                    distInMod += firstIntersection["dist to get to x1y1"]


        rayTracks[rayNum][0].append(["distInFuel",distInFuel])
        rayTracks[rayNum][0].append(["distInMod", distInMod] )
        return distInFuel,distInMod

    else:
        tracks = rayTracks[rayNum]
        mat = rayTracks[0][2]["mat"]
        cellID = rayTracks[0][2]["cellID"]
        whereRayIs = tracks[0][1]
        psi = 0.0;
        if not whereRayIs: psi = sim.q[cellID][0]/(4.0*pi*mat.SigmaT)
        else             : psi = sim.q[cellID][1]/(4.0*pi*mat.SigmaT)
        for i in range(len(tracks)-1):
            mat    = tracks[i+1]["mat"   ]; idVal  = tracks[i+1]["idVal"]
            length = tracks[i+1]["length"]; l      = tracks[i+1]["l"]
            cellID = tracks[i+1]["cellID"]; 
            deltaPsi = (psi - (sim.q[cellID][idVal]/(4.0*pi*mat.SigmaT)))*(1.0-exp(-mat.SigmaT*length))
            psi  -= deltaPsi
            if ( l < 250.0 ):
                sim.phi[cellID][idVal] += 4.0*pi*deltaPsi
        distInFuel = tracks[0][2][1]
        distInMod  = tracks[0][3][1]
        return distInFuel,distInMod
 
    

# Initialize Geometry

#radii = [0.5,0.4,0.3,0.2,0.1]
#radii = [0.4,0.3,0.2,0.1]
#radii = [0.5,0.3]
radii = [0.39128]

cellColors, circleColors = findColors(radii)


# Initialize Materials
#        Sigma T, nuSigma F, Sigma S
mod  = material(3.60, 0.00, 2.10)
fuel = material(0.93, 1.98, 0.30)

xPlanes = [xPlane(-0.63,'ref'),xPlane(0.63,'ref'),xPlane(1.89,'ref'),xPlane(3.15,'ref')]
yPlanes = [yPlane(-0.63,'ref'),yPlane(0.63,'ref'),yPlane(1.89,'ref'),yPlane(3.15,'ref')]

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

allRegionsCell1 = [cells[0]]+cells[0].C
numRegionsPerCell = 1 + len(radii)


volumes = getVolumes(numRegionsPerCell,cells[0].C,sideLength)

random.seed(5)
for i in range(451):
    x = random.random()
    x = random.random()
    x = random.random()
    x = random.random()
    x = random.random()

k_guess = 1.0

q_guess = []
phi_guess = []
oldFissionSource = []
newFissionSource = []
for i in range(len(cells)):
    q_guess.append([0.0]*numRegionsPerCell)
    phi_guess.append([0.0]*numRegionsPerCell)
    oldFissionSource.append([0.0]*numRegionsPerCell)
    newFissionSource.append([0.0]*numRegionsPerCell)

for c in range(len(cells)):
    for i in range(numRegionsPerCell):
        mat = allRegionsCell1[i].mat
        q_guess[c][i] = ( mat.SigmaS + mat.SigmaF )
        oldFissionSource[c][i] = mat.SigmaF

converged = False
firstIteration = True
numRaysPerRun = 20

counter = 0
kVals = []
while not converged:
    #random.seed(5)
    sim = simulation(phi_guess,q_guess,k_guess) 
    if firstIteration:
        distInMod,distInFuel = runRays(sim,numRaysPerRun,firstIteration)
    else:
        runRays(sim,numRaysPerRun,firstIteration)


    #inv_trackLength = 1.0/(250.0*numRaysPerRun)

    inv_trackLength = 1.0/(distInMod+distInFuel)

    for cell in cells:
        sim.phi[cell.id] = [inv_trackLength * x for x in sim.phi[cell.id]]


    # Update phi
    for cell in cells:
        for i in range(numRegionsPerCell):
            mat = allRegionsCell1[i].mat
            idVal = allRegionsCell1[i].id
            sim.phi[cell.id][i] = sim.phi[cell.id][i]/(mat.SigmaT*volumes[idVal]) + sim.q[cell.id][i]/mat.SigmaT
            sim.q[cell.id][i] = ( mat.SigmaS*sim.phi[cell.id][i] + mat.SigmaF*sim.phi[cell.id][i] )
            newFissionSource[cell.id][i] = mat.SigmaF*sim.phi[cell.id][i]

    kNumer = 0
    kDenom = 0
    for i in range(len(newFissionSource)):
        kNumer += sum(newFissionSource[i])
        kDenom += sum(oldFissionSource[i])
    sim.k = kNumer / kDenom

    for cell in cells:
        for i in range(numRegionsPerCell):
            oldFissionSource[cell.id][i] = newFissionSource[cell.id][i]/sim.k
            sim.q[cell.id][i] /= sim.k 



    #sim.k = sum(sum(newFissionSource))/sum(sum(oldFissionSource))

    print()
    print(sim.k)
    #print(" ".join('%0.5f' % item for item in newFissionSource))

    k_guess = sim.k
    kVals.append(sim.k)
    phi_guess = []
    for i in range(len(cells)):
        phi_guess.append([0.0]*numRegionsPerCell)
    q_guess = sim.q[:]
    # Check if converged
    counter += 1
    if counter > 1000: converged = True
    firstIteration = False

#print(rayTracks)
#for cell in cells:
#    plotBox(ax,cell)

#plt.xlim(-1,3.5)
#plt.ylim(-1,3.5)
#plt.show()
plt.plot(kVals)
plt.show()
