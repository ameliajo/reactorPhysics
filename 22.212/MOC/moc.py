import matplotlib.pyplot as plt
import random
from geom import *
from plotting import *
from materials import *

random.seed(3)
fig, ax = plt.subplots() 

class ray:
    def __init__(self,x,y,sin,cos,length):
        self.x0 = x
        self.y0 = y
        self.sin = sin
        self.cos = cos
        self.l = length

def findAllIntersections(xPlanes,yPlanes,circles,r):
    intersections = []

    for xPlane in xPlanes:
        t = (xPlane.x - r.x0)/r.cos
        if ( t > 0 ):
            x_int = round(r.x0+t*r.cos,10)
            y_int = round(r.y0+t*r.sin,10)
            intersections.append({"x-int":x_int,"y-int":y_int,"t":t,"bc":xPlane.bc,"dir":"x"})

    for yPlane in yPlanes:
        t = (yPlane.y - r.y0)/r.sin
        if ( t > 0 ):
            x_int = round(r.x0+t*r.cos,10)
            y_int = round(r.y0+t*r.sin,10)
            intersections.append({"x-int":x_int,"y-int":y_int,"t":t,"bc":yPlane.bc,"dir":"y"})

    for circle in circles:
        circ_x0 = circle.x
        circ_y0 = circle.y
        R = circle.r
        A = r.sin*r.sin+r.cos*r.cos
        B = 2.0*r.x0*r.cos - 2.0*circ_x0*r.cos + 2.0*r.y0*r.sin - 2.0*circ_y0*r.sin
        C = r.x0*r.x0 + circ_x0*circ_x0 + r.y0*r.y0 + circ_y0*circ_y0 - R*R

        plus_minus = [1.0,-1.0]
        for pm in plus_minus:
            t = (-B + pm*(B*B-4*A*C)**0.5)/(2.0*A)
            if (abs(t.imag) < 1.0e-20 and t.real > 0.0):
                t = t.real
                x_int = r.x0 + t*r.cos
                y_int = r.y0 + t*r.sin
                intersections.append({"x-int":x_int,"y-int":y_int,"t":t,"bc":"None","dir":"None"})
    return intersections


def findFirstIntersection(intersections,r):
    t_min = intersections[0]["t"]
    first_intersection = intersections[0]
    for intersection in intersections:
        t = intersection["t"]
        if ( t < t_min ):
            first_intersection = intersection
            t_min = t

    r.x1 = first_intersection["x-int"]
    r.y1 = first_intersection["y-int"]
    dist_traveled = ((r.x1-r.x0)*(r.x1-r.x0)+(r.y1-r.y0)*(r.y1-r.y0))**0.5
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
    currentCircle = potential_circles[0]
    for circle in potential_circles:
        if (circle.r < currentCircle.r):
            currentCircle = circle
    return currentCircle


def isRayInCircle(r,circle):
    smudged_x = r.x0 + 0.001*r.cos
    smudged_y = r.y0 + 0.001*r.sin
    #return (r.x0-circle.x)**2+(r.y0-circle.y)**2 <= (circle.r)**2
    return (smudged_x-circle.x)**2+(smudged_y-circle.y)**2 <= (circle.r)**2


def runRays():
    # Initialize Ray 
    x0 = (cell.R.x-cell.L.x)*(random.random()-0.5)
    y0 = (cell.U.y-cell.D.y)*(random.random()-0.5)
    cos = 2.0*random.random()-1.0
    sin = 2.0*random.random()-1.0

    r = ray(x0,y0,sin,cos,4) # should be 4 length
    
    xPlanes = [cell.L,cell.R]
    yPlanes = [cell.U,cell.D]
    circles = cell.C

    
    counter = 0
    
    while(r.l > 0):
        r.xMax = r.x0+r.cos*r.l
        r.yMax = r.y0+r.sin*r.l
    
        #plotRay(r,colors[counter])
        intersections = findAllIntersections(xPlanes,yPlanes,circles,r)
        if (len(intersections)==0): break
        firstIntersection = findFirstIntersection(intersections,r)
        f = firstIntersection
        #plotRaySegment(r,colors[counter%len(colors)],firstIntersection)
        currentCircle = whichCircleIsRayIn(r,circles)
        if (currentCircle):
            plotRaySegment(r,currentCircle.color,firstIntersection)
        else:
            plotRaySegment(r,'#1e4877',firstIntersection)
        if (firstIntersection["t"] > r.l): 
            r.x = r.xMax
            r.y = r.yMax
            break
        updateRay(r,firstIntersection,counter)
       
        counter += 1
    

# Initialize Geometry
radii = [0.2,0.3,0.4,0.5]
color = ['#ee3e32','#f68838','#fbb021','#1b8a5a']
# box(box_x,box_y,sideLength,radii vec, mod mat, fuel mat, color vec)
cell = box(0.0,0.0,1.0,radii,1,2,color)



for i in range(20):
    runRays()


plotBox(ax,cell)

plt.xlim(-1,1)
plt.ylim(-1,1)
plt.show()

