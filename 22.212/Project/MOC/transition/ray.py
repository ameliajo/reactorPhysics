import numpy as np
import cmath
from numpy.random import random_sample as rand
from math import pi

from region import *
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
        self.x = self.r[0]
        self.y = self.r[1]
        self.ux = self.u[0]
        self.uy = self.u[1]


class Segment():
    def __init__(self, r0, r1, mu, region_num, active=True):
        self.r0 = r0
        self.r1 = r1
        self.mu = mu
        self.region = region_num
        self.d = np.linalg.norm(r1-r0)/mu
        self.active = active


def drawRay(ray, surfaces, regions, deadzone):
    while ray.length < ray.maxDist:

        bestInt = {"t":1e5}

        allSurfaceCrossings = [crossCircle(ray.r,ray.u,circle) for circle in surfaces[0]] +  \
                              [crossXPlane(ray.r,ray.u,xPlane) for xPlane in surfaces[1]] +  \
                              [crossYPlane(ray.r,ray.u,yPlane) for yPlane in surfaces[2]] 

        for surface in allSurfaceCrossings:  # make bestInt equal to current
            for intersection in surface:     # intersection if currents got 
                                             # shorter distance 
                bestInt = intersection if intersection['t'] < bestInt['t'] else bestInt


        r = np.array([bestInt['x'],bestInt['y']])

        fullDistTraveled = ((r[0]-ray.r[0])**2 + (r[1]-ray.r[1])**2)**0.5 / ray.mu
        ray.length += fullDistTraveled


        region_id = 0
        for region in regions:
            if region.evaluate(ray.r):
                region_id = region.uid
                break
 
        if ray.length < deadzone:
            segment = Segment(ray.r, r, ray.mu, region_id, active=False)
        else:
            regions[region_id].activeDist += fullDistTraveled
            ray.active_length += fullDistTraveled
            segment = Segment(ray.r, r, ray.mu, region_id, active=True)

        ray.segments.append(segment)
        ray.r = r

        if bestInt['surface'].boundary_type == 'reflection':
            ray.u = np.array([-ray.u[0], ray.u[1]]) if bestInt['surface'].type == 'x' \
               else np.array([ray.u[0], -ray.u[1]])


        # Move ray forward a small bit to insure location in new region
        smudge = 1e-11
        ray.r += ray.u*smudge

    return ray



