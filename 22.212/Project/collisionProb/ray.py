import numpy as np
import cmath
from numpy.random import random_sample as rand
from math import pi
from crossingGeom import *
from geometry import *

class Neutron():
    def __init__(self,sideLen=0,reg=None,init=True):
        if init and not reg:
            self.r = np.array([rand(),rand()])*sideLen
        elif init and reg.type == 'circle':
            l, r = reg.x0 - reg.r, reg.x0 + reg.r
            d, u = reg.y0 - reg.r, reg.y0 + reg.r
            inCircle = False
            while not inCircle:
               self.r = np.array([l+rand()*(r-l),d+rand()*(u-d)])
               inCircle = (self.r[0]-reg.x0)**2 + (self.r[1]-reg.y0)**2 < reg.r**2

        if init:
            self.mu = np.cos((2.0*rand()-1.0)*pi*0.5)
            theta = rand()*2*pi
            self.u = np.array([np.cos(theta), np.sin(theta)])


            


class Segment():
    def __init__(self, ray, nextStop, regionID, active):
        self.start = ray.r
        self.end   = nextStop
        self.mu = ray.mu
        self.regionID = regionID 
        self.d = np.linalg.norm(self.start-self.end)/self.mu
        self.active = active


def advance(n, surfaces, regions):

    bestInt = {"t":1e5}

    allSurfaceCrossings = [crossCircle(n,circle) for circle in surfaces[0]] +  \
                              [crossXPlane(n,xPlane) for xPlane in surfaces[1]] +  \
                              [crossYPlane(n,yPlane) for yPlane in surfaces[2]] 

    for surface in allSurfaceCrossings:  # make bestInt equal to current
        for intersection in surface:     # intersection if currents got 
                                             # shorter distance 
            bestInt = intersection if intersection['t'] < bestInt['t'] else bestInt


    r = np.array([bestInt['x'],bestInt['y']])

    fullDistTraveled = ((r[0]-n.r[0])**2 + (r[1]-n.r[1])**2)**0.5 / n.mu

    regionID = [region for region in regions if region.evaluate(n.r)][0].ID

    #n.r = r
    n_next = Neutron()
    n_next.r = r

        
    surfaceWeHit = bestInt['surface']
    if surfaceWeHit.BC == 'reflection':
        n_next.u = np.array([-n.u[0], n.u[1]]) if isinstance(surfaceWeHit,XPlane) \
              else np.array([n.u[0], -n.u[1]])



    n_next.r += n.u*1e-11

    return n_next



