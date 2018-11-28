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

    def evaluate(cls, r):
        """Check if r is located within the region
        r : 2-tuple coordinates of point
        """
        for idx, surface in enumerate(cls.surfaces):
            sgn = cls.orientations[idx]
            if sgn*surface.evaluate(r) < 0:
                return False
        return True

        
    
