import numpy as np

class Surface(object):
    def __init__(self, surface_id, boundary_type, name =''):
        self.id = surface_id
        self.boundary_type = boundary_type
        self.name = name

class XPlane():
    type = 'x'
    def __init__(self, surface_id, x0, boundary_type='transmission',name=''):
        self.x0 = x0
        self.id = surface_id
        self.boundary_type = boundary_type
        self.name = name



    def evaluate(self, point):
        return point[0] - self.x0


class YPlane():
    type = 'y'
    def __init__(self, surface_id, y0, boundary_type='transmission',name=''):
        self.y0 = y0
        self.id = surface_id
        self.boundary_type = boundary_type
        self.name = name



    def evaluate(self, point):
        return point[1] - self.y0

class Circle():
    type = 'circle'
    def __init__(self, surface_id, x0, y0, R, name='',boundary_type='transmission'):
        self.x0 = x0
        self.y0 = y0
        self.r = R
        self.id = surface_id
        self.boundary_type = boundary_type
        self.name = name



    def evaluate(self, point):
        x = point[0] - self.x0
        y = point[1] - self.y0
        return x**2 + y**2 - self.r**2

    def rayIsInCircle(self,r):
        x = r[0] - self.x0
        y = r[1] - self.y0
        return x**2 + y**2 < self.r**2





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



        
 

    
