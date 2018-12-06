import numpy as np
from runMC     import *
from geometry  import *
from materials import *

pitch    = 1.26
sideLen  = pitch*3
radius   = 0.39218
numRays  = 100

circX = [pitch*(i+0.5) for i in range(3)]
circY = [pitch*(i+0.5) for i in range(3)]

boundaries = ['reflection', 'transmission', 'transmission', 'reflection']
XPlanes = [XPlane(i*pitch, boundaries[i]) for i in range(4)]
YPlanes = [YPlane(i*pitch, boundaries[i]) for i in range(4)]


circles = [ Circle( circX[i], circY[j], radius ) for j in range(3) for i in range(3) ]
fuelVec = [ Region( [ (circles[i], False) ], fuelClass, 'fuel' ) for i in range(9) ] 
modVec  = [ Region( [ (XPlanes[i], True), (XPlanes[i+1], False),               \
                      (YPlanes[j], True), (YPlanes[j+1], False),               \
                      (circles[-1], True) ], modClass, 'mod')              \
                      for j in range(3) for i in range(3) ]

regions  = fuelVec + modVec

for i, region in enumerate(regions):
    region.ID = i

surfaces = [circles,XPlanes,YPlanes]
runMC(numRays, surfaces, regions, sideLen, False)


