import numpy as np

class Surface(object):
    def __init__(self, surface_id, boundary_type, name =''):
        self.id = surface_id
        self.boundary_type = boundary_type
        self.name = name

class XPlane(Surface):
    type = 'x-plane'

    def __init__(self, surface_id, x0, boundary_type='transmission',name=''):
        super().__init__(surface_id, boundary_type, name=name)
        self.x0 = x0

    def evaluate(self, point):
        return point[0] - self.x0

class YPlane(Surface):
    type = 'y-plane'

    def __init__(self, surface_id, y0, boundary_type='transmission',name=''):
        super().__init__(surface_id, boundary_type, name=name)
        self.y0 = y0

    def evaluate(self, point):
        return point[1] - self.y0

class Circle(Surface):
    type = 'circle'

    def __init__(self, surface_id, x0, y0, R, name='',boundary_type='transmission'):
        super().__init__(surface_id, boundary_type, name=name)
        self.x0 = x0
        self.y0 = y0
        self.r = R

    def evaluate(self, point):
        x = point[0] - self.x0
        y = point[1] - self.y0
        return x**2 + y**2 - self.r**2






    
