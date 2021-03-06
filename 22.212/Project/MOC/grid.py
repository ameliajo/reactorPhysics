import numpy as np
import sys
from runRays   import *
from geometry  import *
from materials import *

verbose = True
sphIter = 1
if len(sys.argv) > 1:
    if sys.argv[1] == 'quiet': verbose = False
if len(sys.argv) > 2:
    sphIter = int(sys.argv[2])
        

pitch    = 1.26
sideLen  = pitch*3
radius   = 0.39218
numRays  = 100
rayLen   = 100.0
deadzone = 10

circX = [pitch*(i+0.5) for i in range(3)]
circY = [pitch*(i+0.5) for i in range(3)]

boundaries = ['reflection', 'transmission', 'transmission', 'reflection']
XPlanes = [XPlane(i*pitch, boundaries[i]) for i in range(4)]
YPlanes = [YPlane(i*pitch, boundaries[i]) for i in range(4)]

ngroup = 9
phiInitF = np.ones([ngroup,])
phiInitM = np.ones([ngroup,])*0.1

circles = [ Circle( circX[i], circY[j], radius ) for j in range(3) for i in range(3) ]
fuelVec = [ Region( [ (circles[i], False) ], [fuelClassHi,fuelClassLo][i%2], phiInitF ) for i in range(9) ] 
modVec  = [ Region( [ (XPlanes[i], True), (XPlanes[i+1], False),               \
                      (YPlanes[j], True), (YPlanes[j+1], False),               \
                      (circles[-1], True) ], modClass, phiInitM )              \
                      for j in range(3) for i in range(3) ]

regions  = fuelVec + modVec

for i, region in enumerate(regions):
    region.ID = i

surfaces = [circles,XPlanes,YPlanes]
runRays(numRays, surfaces, regions, sideLen, ngroup, False, rayLen, deadzone, verbose, sphIter)


