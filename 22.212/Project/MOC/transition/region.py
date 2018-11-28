import numpy as np

class Region():
    def __init__(self, surfaces, uid, mat,phi=0, orientations = [1, -1, -1, 1, 1],vol=0):
        self.surfaces = surfaces
        self.orientations = orientations
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
            sgn = self.orientations[idx]
            if sgn*surface.evaluate(r) < 0:
                return False
        return True


    def isInRegion(self, r):
        for surface in self.surfaces:
            if surface.type == 'circle':
                return surface.rayIsInCircle(r)



        
    
