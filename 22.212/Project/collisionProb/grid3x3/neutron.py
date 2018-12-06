import numpy as np
import cmath
from numpy.random import random_sample as rand
from math import pi
from geometry import *


def crossCircle(ray, circle):
    c = circle
    x,   y = ray.r
    ux, uy = ray.u
            
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

def crossXPlane(ray,xPlane):
    x,   y = ray.r
    ux, uy = ray.u
    t = (xPlane.x0-x)/ux
    return {"x":x+t*ux,"y":y+t*uy,"t":t,"surface":xPlane} if t > 0 else None


def crossYPlane(ray,yPlane):
    x,   y = ray.r
    ux, uy = ray.u
    t = (yPlane.y0-y)/uy
    return {"x":x+t*ux,"y":y+t*uy,"t":t,"surface":yPlane} if t > 0 else None


class Neutron():
    def __init__(self,sideLen=0,reg=None):
        if not reg:
            self.r = np.array([rand(),rand()])*sideLen
        elif reg.type == 'circle':
            l, r = reg.x0 - reg.r, reg.x0 + reg.r
            d, u = reg.y0 - reg.r, reg.y0 + reg.r
            inCircle = False
            while not inCircle:
               self.r = np.array([l+rand()*(r-l),d+rand()*(u-d)])
               inCircle = (self.r[0]-reg.x0)**2 + (self.r[1]-reg.y0)**2 < reg.r**2

        self.mu = np.cos((2.0*rand()-1.0)*pi*0.5)
        theta = rand()*2*pi
        self.u = np.array([np.cos(theta), np.sin(theta)])



def weCollide(distTraveled,mat):
    randomSampleDist = -np.log(rand())/mat.SigT
    return randomSampleDist < distTraveled



def advance(n, surfaces, regions):

    bestInt = {"t":1e5}

    allSurfaceCrossings = [crossCircle(n,circle) for circle in surfaces[0]] +  \
                          [crossXPlane(n,xPlane) for xPlane in surfaces[1]] +  \
                          [crossYPlane(n,yPlane) for yPlane in surfaces[2]] 

    for intersection in allSurfaceCrossings:
        if intersection and intersection['t'] < bestInt['t']:
            bestInt = intersection

    r = np.array([bestInt['x'],bestInt['y']])

    fullDistTraveled = ((r[0]-n.r[0])**2 + (r[1]-n.r[1])**2)**0.5 / n.mu

    for region in regions:
        if region.evaluate(n.r): 
            regionID = region.ID
            break

    surfaceWeHit = bestInt['surface']
    if surfaceWeHit.BC == 'reflection':
        n.u = np.array([-n.u[0], n.u[1]]) if isinstance(surfaceWeHit,XPlane) \
         else np.array([n.u[0], -n.u[1]])

    n.r = r + n.u*1e-11

    return fullDistTraveled,regionID



