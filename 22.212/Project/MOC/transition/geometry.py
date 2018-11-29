import numpy as np

class Surface(object):
    def __init__(self, surface_id, BC):
        self.id = surface_id
        self.BC = BC

class XPlane():
    type = 'x'
    def __init__(self, surface_id, x0, BC='transmission'):
        self.x0 = x0
        self.id = surface_id
        self.BC = BC


    def areWeInHere(xPlane, x, y, onPosSide):
        return (x > xPlane.x0) == onPosSide


class YPlane():
    type = 'y'
    def __init__(self, surface_id, y0, BC='transmission'):
        self.y0 = y0
        self.id = surface_id
        self.BC = BC

    def areWeInHere(yPlane, x, y, onPosSide):
        return (y > yPlane.y0) == onPosSide

class Circle():
    type = 'circle'
    def __init__(self, surface_id, x0, y0, R, BC='transmission'):
        self.x0 = x0
        self.y0 = y0
        self.r = R
        self.id = surface_id
        self.BC = BC

    def areWeInHere(circle, x, y, onPosSide):
        return ((x-circle.x0)**2 + (y-circle.y0)**2 > circle.r**2) == onPosSide


class Region():
    def __init__(self, surfacesAndSides, uid, mat, phi=0, vol=0):
        surfacesTHENsides= list(zip(*surfacesAndSides))
        self.surfaces     = surfacesTHENsides[0]
        self.orientations = surfacesTHENsides[1]
        self.uid = uid
        self.phi = phi
        self.vol = vol
        self.mat = mat
        nGroups = len(phi)

        self.q = np.zeros([nGroups,])
        self.tracks_phi = np.zeros([nGroups,])
        self.activeDist = 0.0

    def evaluate(self, r):
        for idx, surface in enumerate(self.surfaces):
            onPosSide = self.orientations[idx]
            if not surface.areWeInHere(r[0],r[1],onPosSide):
                return False

        return True




 

    
