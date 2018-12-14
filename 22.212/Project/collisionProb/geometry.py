
class Material:
    def __init__(self,SigT,name):
        self.SigT = SigT
        self.name = name



def makeGeometry(pitch,radius,hole):
    circX = [pitch*(i+0.5) for i in range(3)]
    circY = [pitch*(i+0.5) for i in range(3)]

    XPlanes = [XPlane(i*pitch) for i in range(4)]
    YPlanes = [YPlane(i*pitch) for i in range(4)]

    circles = [Circle(circX[i],circY[j],radius) for j in range(3) for i in range(3)]
    fRegs = [Region([ (circles[i],False) ],['fuel (high enr.)','fuel (low enr.)'][i%2]) for i in range(9)] 

    if hole: circles.pop(4)

    mRegs = [Region([ (XPlanes[i],True), (XPlanes[i+1],False),        \
                      (YPlanes[j],True), (YPlanes[j+1],False),        \
                      (circles[-1],True) ],                           \
                      'mod')                                     \
                     for j in range(3) for i in range(3) ]

    regions  = fRegs + mRegs
    
    for i, region in enumerate(regions): region.ID = i
    surfaces = [circles,XPlanes,YPlanes]

    return regions,surfaces


class XPlane():
    type = 'x'
    def __init__(self, x0):
        self.x0 = x0
    def areWeInHere(xPlane, x, y, onPosSide):
        return (x > xPlane.x0) == onPosSide

class YPlane():
    type = 'y'
    def __init__(self, y0):
        self.y0 = y0
    def areWeInHere(yPlane, x, y, onPosSide):
        return (y > yPlane.y0) == onPosSide

class Circle():
    type = 'circle'
    def __init__(self, x0, y0, R):
        self.x0, self.y0 = x0, y0
        self.r = R
    def areWeInHere(circle, x, y, onPosSide):
        return ((x-circle.x0)**2 + (y-circle.y0)**2 > circle.r**2) == onPosSide

class Region():
    def __init__(self, surfacesAndSides, name):
        self.name = name
        surfacesTHENsides = list(zip(*surfacesAndSides))
        self.surfaces     = surfacesTHENsides[0]
        self.orientations = surfacesTHENsides[1]

    def evaluate(self, x,y):
        for i, surface in enumerate(self.surfaces):
            onPosSide = self.orientations[i]
            if not surface.areWeInHere(x,y,onPosSide): return False
        return True




 

    
