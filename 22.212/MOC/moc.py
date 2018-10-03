import matplotlib.pyplot as plt
import random
from geom import *
import matplotlib.cm as cm
import numpy as np

random.seed(3)
fig, ax = plt.subplots() 

class ray:
    def __init__(self,x,y,sin,cos,length):
        self.x0 = x
        self.y0 = y
        self.sin = sin
        self.cos = cos
        self.l = length

def plotIntersection(r,t):
    x_int = r.x0 + t * r.cos
    y_int = r.y0 + t * r.sin
    #plt.plot(x_int,y_int,'bo')
    return (x_int,y_int)

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






def findAllIntersections(xPlanes,yPlanes,r):
    intersections = []
    for xPlane in xPlanes:
        t = (xPlane.x - r.x0)/r.cos
        if ( t > 0 ):
            x_int,y_int = plotIntersection(r,t)
            intersections.append({"x-int":x_int,"y-int":y_int,"t":t,"bc":xPlane.bc,"dir":"x"})

    for yPlane in yPlanes:
        t = (yPlane.y - r.y0)/r.sin
        if ( t > 0 ):
            x_int,y_int = plotIntersection(r,t)
            intersections.append({"x-int":x_int,"y-int":y_int,"t":t,"bc":yPlane.bc,"dir":"y"})
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

def updateRay(r,firstIntersection):
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





circle1 = circle(0,0,0.5,'U')
box1 = box(xPlane(-1,'ref'),xPlane(1,'ref'),
           yPlane(1,'ref'),yPlane(-1,'ref'),'H',circle1)

x0 = (box1.R.x-box1.L.x)*(random.random()-0.5)
y0 = (box1.U.y-box1.D.y)*(random.random()-0.5)
cos = 2.0*random.random()-1.0
sin = 2.0*random.random()-1.0
r = ray(x0,y0,sin,cos,50) # should be 4 length


#plotRay(r,"r")

xPlanes = [box1.L,box1.R]
yPlanes = [box1.U,box1.D]


counter = 0
colors = ["crimson","deeppink","fuchsia","violet","darkviolet","rebeccapurple","slateblue","slategrey","dodgerblue","skyblue","cadetblue","c","darkslategray","mediumturquoise","mediumaquamarine","mediumseagreen","g","lightgreen","lawngreen", "yellowgreen", "olive", "lemonchiffon", "darkgoldenrod", "orange", "navajowhite", "darkorange", "peachpuff", "coral", "mistyrose", "maroon", "lightcoral", "lightgray", "grey", "k"]
while(r.l > 0):
    r.xMax = r.x0+r.cos*r.l
    r.yMax = r.y0+r.sin*r.l

    #plotRay(r,colors[counter])
    intersections = findAllIntersections(xPlanes,yPlanes,r)
    if (len(intersections)==0): break
    firstIntersection = findFirstIntersection(intersections,r)
    plotRaySegment(r,colors[counter],firstIntersection)
    if (firstIntersection["t"] > r.l): 
        r.x = r.xMax
        r.y = r.yMax
        break
    updateRay(r,firstIntersection)
    
    counter += 1


#intersections = findAllIntersections(xPlanes,yPlanes,r)
#firstIntersection = findFirstIntersection(intersections)
#plt.plot(firstIntersection["x-int"],firstIntersection["y-int"],'go')
#updateRay(r,firstIntersection)
#plotRay(r,"g")


box1.plot(ax)





plt.xlim(-4,4)
plt.ylim(-4,4)
plt.show()



