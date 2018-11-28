import numpy as np
from main import main 

from surface import XPlane, YPlane, Circle
from tools import *
from plotting import *
from region import Region
from ray import * 
from mat import *

pitch = 1.26
sideLen= 1.26*3
radius = 0.39218
n_rays = 100
rayLen = 100.0

# Surfaces starting with 100 will be circles
# Surfaces starting with 200 will be planes
circX = [pitch*(i+0.5) for i in range(3)]
circY = [pitch*(i+0.5) for i in range(3)]


boundaries = ['reflection', 'transmission', 'transmission', 'reflection']
XPlanes = [XPlane(201+i, i*pitch,     boundaries[i]) for i in range(4)]
YPlanes = [YPlane(205+i, (3-i)*pitch, boundaries[i]) for i in range(4)]


# Regions
ngroup = 10
fuel_phi_guess = np.ones([ngroup,])
mod_phi_guess = np.ones([ngroup,])*0.1

circles = []
modVec = []
counter = 0
for j in range(3):
    for i in range(3):
        circles.append(Circle(101+counter, circX[i], circY[j], radius))
        modVec.append(Region([XPlanes[i], XPlanes[i+1], YPlanes[j], YPlanes[j+1], 
                      circles[-1]],counter+9, modClass, mod_phi_guess))

        counter += 1



regions = [Region([circles[i]], i, fuelClass,fuel_phi_guess,[-1]) for i in range(9)]
regions += modVec
surfaces = [circles,XPlanes,YPlanes]
#main(n_rays, surfaces, regions, sideLen, ngroup, False, rayLen)
main(n_rays, surfaces, regions, sideLen, ngroup, True, rayLen)

