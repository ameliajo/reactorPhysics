import numpy as np
from main      import *
from geometry  import *
from materials import *

pitch   = 1.26
sideLen = pitch*3
radius  = 0.39218
numRays = 100
rayLen  = 100.0

circX = [pitch*(i+0.5) for i in range(3)]
circY = [pitch*(i+0.5) for i in range(3)]

boundaries = ['reflection', 'transmission', 'transmission', 'reflection']
XPlanes = [XPlane(100+i, i*pitch,     boundaries[i]) for i in range(4)]
YPlanes = [YPlane(104+i, (3-i)*pitch, boundaries[i]) for i in range(4)]

ngroup = 10
phiInitF = np.ones([ngroup,])
phiInitM = np.ones([ngroup,])*0.1

circles, modVec = [], []
count = 0
for j in range(3):
    for i in range(3):
        circles.append(Circle(count, circX[i], circY[j], radius))
        modVec.append(Region([(XPlanes[i],True ),(XPlanes[i+1],False), 
                              (YPlanes[j],False),(YPlanes[j+1],True ), 
                              (circles[-1],True)],count+9, modClass, phiInitM))
        count += 1


regions = [Region([(circles[i],False)], i, fuelClass,phiInitF) for i in range(9)]
regions += modVec
surfaces = [circles,XPlanes,YPlanes]
runRays(numRays, surfaces, regions, sideLen, ngroup, False, rayLen)


