import numpy as np
import cmath
from numpy.random import random_sample as rand
from math import pi
from geometry import *


def crossCircle(x,y,ux,uy, circle):
    c = circle
            
    A = ux**2 + uy**2
    B = 2.0*(x*ux - c.x0*ux + y*uy - c.y0*uy)
    C = x**2 - 2.0*x*c.x0 + c.x0*c.x0 + y**2 - 2.0*y*c.y0 + c.y0**2 - c.r**2

    intersections = None
    best_t = 1e5
    for pm in [1.0,-1.0]:
        if B*B - 4*A*C < 0: continue
        t = (-B + pm*(B*B-4*A*C)**0.5)/(2*A)
        if (abs(t.imag) < 1.0e-20 and t.real > 0.0 and t.real < best_t):
            t = t.real
            best_t = t
            xInt, yInt = x + t*ux, y + t*uy
            intersections = {"x":xInt,"y":yInt,"t":t,"surface":circle}

    return intersections

def crossXPlane(x,y,ux,uy,xPlane):
    t = (xPlane.x0-x)/ux
    return {"x":x+t*ux,"y":y+t*uy,"t":t,"surface":xPlane} if t > 0 else None


def crossYPlane(x,y,ux,uy,yPlane):
    t = (yPlane.y0-y)/uy
    return {"x":x+t*ux,"y":y+t*uy,"t":t,"surface":yPlane} if t > 0 else None




def weCollide(distTraveled,SigT):
    return -np.log(rand())/SigT < distTraveled



def advance(x,y,ux,uy,mu, surfaces, regions):
    # Don't need to differentiate between different moderator regions
    allSurfaceCrossings = [crossCircle(x,y,ux,uy,circle) for circle in  surfaces[0]] +  \
                          [crossXPlane(x,y,ux,uy,xPlane) for xPlane in [surfaces[1][0],\
                                                                surfaces[1][3]]] + \
                          [crossYPlane(x,y,ux,uy,yPlane) for yPlane in [surfaces[2][0],\
                                                                surfaces[2][3]]]

    bestInt = {"t":1e5}
    for intersection in allSurfaceCrossings:
        if intersection and intersection['t'] < bestInt['t']:
            bestInt = intersection


    fullDistTraveled = ((bestInt['x']-x)**2 + (bestInt['y']-y)**2)**0.5 / mu

    for region in regions:
        if region.evaluate(x,y): 
            regionID = region.ID
            break

    ux = -ux if bestInt['surface'].type == 'x' else ux
    uy = -uy if bestInt['surface'].type == 'y' else uy

    x = bestInt['x'] + ux*1e-11
    y = bestInt['y'] + uy*1e-11

    return fullDistTraveled,regionID,x,y,ux,uy



