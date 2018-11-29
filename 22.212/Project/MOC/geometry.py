import numpy as np

class XPlane():
    def __init__(self, x0, BC='transmission'):
        self.x0 = x0
        self.BC = BC

    def areWeInHere(xPlane, x, y, onPosSide):
        return (x > xPlane.x0) == onPosSide



class YPlane():
    type = 'y'
    def __init__(self, y0, BC='transmission'):
        self.y0 = y0
        self.BC = BC

    def areWeInHere(yPlane, x, y, onPosSide):
        return (y > yPlane.y0) == onPosSide



class Circle():
    type = 'circle'
    def __init__(self, x0, y0, R, BC='transmission'):
        self.x0 = x0
        self.y0 = y0
        self.r = R
        self.BC = BC

    def areWeInHere(circle, x, y, onPosSide):
        return ((x-circle.x0)**2 + (y-circle.y0)**2 > circle.r**2) == onPosSide



class Region():
    def __init__(self, surfacesAndSides, mat, phi):
        surfacesTHENsides = list(zip(*surfacesAndSides))
        self.surfaces     = surfacesTHENsides[0]
        self.orientations = surfacesTHENsides[1]
        self.phi = phi
        self.mat = mat
        self.q = np.zeros([len(phi)])
        self.activeDist = 0.0

    def evaluate(self, r):
        for i, surface in enumerate(self.surfaces):
            onPosSide = self.orientations[i]
            if not surface.areWeInHere(r[0],r[1],onPosSide):
                return False
        return True




 

    
