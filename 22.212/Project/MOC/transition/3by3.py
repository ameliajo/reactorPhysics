import numpy as np
from main import main 

from surface import XPlane, YPlane, Circle
from tools import *
from plotting import *
from region import Region
from ray import Ray, make_segments
from mat import *

pitch = 1.26
length = 1.26*3
radius = 0.39218
# Surfaces starting with 100 will be circles
# Surfaces starting with 200 will be planes
xcoords = [pitch*(i+0.5) for i in range(3)]
ycoords = [pitch*(i+0.5) for i in range(3)]

circles = []
counter = 0
for j in range(3):
    for i in range(3):
        circle = Circle(surface_id=101+counter, x0=xcoords[i], y0=ycoords[j], R=radius)
        counter += 1
        circles.append(circle)


boundaries = ['reflection', 'transmission', 'transmission', 'reflection']
XPlanes = [XPlane(201+i, i*pitch,     boundaries[i]) for i in range(4)]
YPlanes = [YPlane(205+i, (3-i)*pitch, boundaries[i]) for i in range(4)]


# Regions
ngroup = 10
fuel_phi_guess = np.ones([ngroup,])
mod_phi_guess = np.ones([ngroup,])*0.1
regions = []
for i in range(9):
    regions.append(Region([circles[i]], i, fuelClass,'fuel',fuel_phi_guess,[-1]))

modVec = []
counter = 0



for j in range(3):
    for i in range(3):
        surfaces = [XPlanes[i], XPlanes[i+1], YPlanes[j], YPlanes[j+1], circles[counter]]
        modVec.append(Region(surfaces,counter+9, modClass,'mod', mod_phi_guess))
        counter += 1
regions += modVec

surfaces = [circles,XPlanes,YPlanes]

n_rays = 100
main(n_rays, surfaces, regions, length, ngroup, False)

