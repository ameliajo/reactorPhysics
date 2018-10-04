import matplotlib.pyplot as plt
import random
from geom import *

random.seed(1)
#random.seed(2)
#random.seed(3)
#random.seed(4)
#random.seed(5)
#random.seed(6)
#random.seed(7)
fig, ax = plt.subplots() 

class ray:
    def __init__(self,x,y,sin,cos,length):
        self.x0 = x
        self.y0 = y
        self.sin = sin
        self.cos = cos
        self.l = length

def plotRay(r,color):
    x = [r.x0,r.xMax]
    y = [r.y0,r.yMax]
    plt.plot(r.x0,r.y0,color,marker='o')
    plt.plot(x,y,color)
    if (color == "rebeccapurple"):
        plt.plot(r.xMax,r.yMax,"aqua",marker='o')

def plotRaySegment(r,color,firstIntersection):
    x = [r.x0,firstIntersection["x-int"]]
    y = [r.y0,firstIntersection["y-int"]]
    plt.plot(x,y,color)






def findAllIntersections(xPlanes,yPlanes,circles,r):
    intersections = []
    for xPlane in xPlanes:
        t = (xPlane.x - r.x0)/r.cos
        if ( t > 0 ):
            x_int = round(r.x0 + t * r.cos,10)
            y_int = round(r.y0 + t * r.sin,10)
            intersections.append({"x-int":x_int,"y-int":y_int,"t":t,"bc":xPlane.bc,"dir":"x"})

    for yPlane in yPlanes:
        t = (yPlane.y - r.y0)/r.sin
        if ( t > 0 ):
            x_int = round(r.x0 + t * r.cos,10)
            y_int = round(r.y0 + t * r.sin,10)
            intersections.append({"x-int":x_int,"y-int":y_int,"t":t,"bc":yPlane.bc,"dir":"y"})

    for circle in circles:
        circ_x0 = box1.C.x
        circ_y0 = box1.C.y
        R = box1.C.r
        A = r.sin*r.sin+r.cos*r.cos
        B = 2.0*r.x0*r.cos - 2.0*circ_x0*r.cos + 2.0*r.y0*r.sin - 2.0*circ_y0*r.sin
        C = r.x0*r.x0 + circ_x0*circ_x0 + r.y0*r.y0 + circ_y0*circ_y0 - R*R

        t = ( -B + (B*B-4*A*C)**0.5 )/(2.0*A)
        if (abs(t.imag) < 1.0e-20 and t.real > 0.0):
            t = t.real
            x_int = r.x0 + t*r.cos
            y_int = r.y0 + t*r.sin
            intersections.append({"x-int":x_int,"y-int":y_int,"t":t,"bc":"None","dir":"None"})
            #plt.plot(x_int,y_int,"yo")

        t = ( -B - (B*B-4*A*C)**0.5 )/(2.0*A)
        if (abs(t.imag) < 1.0e-20 and t.real > 0.0):
            t = t.real
            x_int = r.x0 + t*r.cos
            y_int = r.y0 + t*r.sin
            intersections.append({"x-int":x_int,"y-int":y_int,"t":t,"bc":"None","dir":"None"})
            #plt.plot(x_int,y_int,"go")

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

def updateRay(r,firstIntersection,counter):
    # do we bounce off ?
    if firstIntersection["bc"] == "ref":
        if   firstIntersection["dir"] == "x":
            r.cos = -r.cos
        elif firstIntersection["dir"] == "y":
            r.sin = -r.sin
        else:
            raise ValueError('Dont know where to reflect to')
   
    dist_traveled = firstIntersection["dist to get to x1y1"]
    # Update position and distance to travel
    r.x0 = r.x1
    r.y0 = r.y1
    r.l = r.l - dist_traveled

    if (r.sin > 0.0): r.y0 += 1e-8
    if (r.cos > 0.0): r.x0 += 1e-8
    if (r.sin <= 0.0): r.y0 -= 1e-8
    if (r.cos <= 0.0): r.x0 -= 1e-8





circle1 = circle(0,0,0.5,'U')
box1 = box(xPlane(-1,'ref'),xPlane(1,'ref'),
           yPlane(1,'ref'),yPlane(-1,'ref'),'H',circle1)

x0 = (box1.R.x-box1.L.x)*(random.random()-0.5)
y0 = (box1.U.y-box1.D.y)*(random.random()-0.5)
cos = 2.0*random.random()-1.0
sin = 2.0*random.random()-1.0
r = ray(x0,y0,sin,cos,20) # should be 4 length

plt.plot(x0,y0,'ro')

xPlanes = [box1.L,box1.R]
yPlanes = [box1.U,box1.D]
circles = [box1.C]


counter = 0

while(r.l > 0):
    r.xMax = r.x0+r.cos*r.l
    r.yMax = r.y0+r.sin*r.l

    #plotRay(r,colors[counter])
    intersections = findAllIntersections(xPlanes,yPlanes,circles,r)
    if (len(intersections)==0): break
    firstIntersection = findFirstIntersection(intersections,r)
    f = firstIntersection
    plotRaySegment(r,colors[counter],firstIntersection)
    if (firstIntersection["t"] > r.l): 
        r.x = r.xMax
        r.y = r.yMax
        break
    updateRay(r,firstIntersection,counter)
    plt.plot(r.x0,r.y0,"bo")
   
    counter += 1


box1.plot(ax)




plt.xlim(-2,2)
plt.ylim(-2,2)
plt.show()

