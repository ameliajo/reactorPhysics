import numpy as np
import cmath

from region import what_region
from crossingGeom import *

def make_segments(ray, surfaces, regions, cutoff_length=300, deadzone=50):
    while ray.length < cutoff_length:

        region_id = what_region(ray.r, regions)
        bestInt = {"t":np.Inf}

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


        if ray.length < deadzone:
            segment = Segment(ray.r, r, ray.mu, region_id, active=False)
        else:
            regions[region_id].tot_track_length += fullDistTraveled
            ray.active_length += fullDistTraveled
            segment = Segment(ray.r, r, ray.mu, region_id, active=True)

        ray.segments.append(segment)
        ray.r = r

        if bestInt['surface'].boundary_type == 'reflection':
            ray.u = np.array([-ray.u[0], ray.u[1]]) if bestInt['surface'].type == 'x-plane' \
               else np.array([ray.u[0], -ray.u[1]])


        # Move ray forward a small bit to insure location in new region
        smudge = 1e-11
        ray.r += ray.u*smudge

    return ray


class Ray():
    def __init__(self, r, theta, varphi):
        self.r = r
        self.u = np.array([np.cos(theta), np.sin(theta)])
        self.varphi = varphi
        self.mu = np.cos(varphi)
        self.region = None
        self.segments = []
        self.length = 0
        self.active_length = 0


class Segment():
    def __init__(self, r0, r1, mu, region_num, active=True):
        self.r0 = r0
        self.r1 = r1
        self.mu = mu
        self.region = region_num

        self.d = np.linalg.norm(r1-r0)/mu
        self.active = active

