import numpy as np
import cmath
from numpy.random import random_sample as rand
from math import pi
from crossingGeom import *

class Ray():
    def __init__(self,sideLen,maxDist):
        self.r = np.array([rand(),rand()])*sideLen
        self.mu = np.cos((2.0*rand()-1.0)*pi*0.5)
        theta = rand()*2*pi
        self.u = np.array([np.cos(theta), np.sin(theta)])
        self.segments = []
        self.length = 0
        self.active_length = 0
        self.maxDist = maxDist 


class Segment():
    def __init__(self, ray, nextStop, region_num, active):
        self.start = ray.r
        self.end   = nextStop
        self.mu = ray.mu
        self.region = region_num
        self.d = np.linalg.norm(self.start-self.end)/self.mu
        self.active = active


def drawRay(ray, surfaces, regions, deadzone):
    while ray.length < ray.maxDist:

        bestInt = {"t":1e5}

        allSurfaceCrossings = [crossCircle(ray,circle) for circle in surfaces[0]] +  \
                              [crossXPlane(ray,xPlane) for xPlane in surfaces[1]] +  \
                              [crossYPlane(ray,yPlane) for yPlane in surfaces[2]] 

        for surface in allSurfaceCrossings:  # make bestInt equal to current
            for intersection in surface:     # intersection if currents got 
                                             # shorter distance 
                bestInt = intersection if intersection['t'] < bestInt['t'] else bestInt


        r = np.array([bestInt['x'],bestInt['y']])

        fullDistTraveled = ((r[0]-ray.r[0])**2 + (r[1]-ray.r[1])**2)**0.5 / ray.mu
        ray.length += fullDistTraveled

        regionID = [region for region in regions if region.evaluate(ray.r)][0].uid
 
        if ray.length < deadzone:
            ray.segments.append(Segment(ray, r, regionID, False))
        else:
            regions[regionID].activeDist += fullDistTraveled
            ray.active_length            += fullDistTraveled
            ray.segments.append(Segment(ray, r, regionID, True))

        ray.r = r

        if bestInt['surface'].BC == 'reflection':
            ray.u = np.array([-ray.u[0], ray.u[1]]) if bestInt['surface'].type == 'x' \
               else np.array([ray.u[0], -ray.u[1]])


        ray.r += ray.u*1e-11

    return ray



