from math import pi
import openmc
import openmc.mgxs
import openmc.model
import numpy as np


radius_fuel = 0.3922
pitch = 1.26

# Basic materials
uo2 = openmc.Material(name='fuel')
uo2.add_element('U', 1, enrichment=3.2)
uo2.add_element('O', 2)
#uo2.add_element('Gd', 0.0007)
uo2.set_density('g/cc', 10.341)


water = openmc.Material(3, "h2o")
water.add_element('H', 2.0)
water.add_element('O', 1.0)
water.set_density('g/cm3', 1.0)

"""
waters = []
for i in range(9):
    thiswater = openmc.Material(30+i, "h2o")
    thiswater.add_element('H', 2.0)
    thiswater.add_element('O', 1.0)
    thiswater.set_density('g/cm3', 1.0)
    waters.append(thiswater)

uo2s = []
for i in range(9):
    thisuo2 = openmc.Material(name='fuel')
    thisuo2.add_element('U', 1, enrichment=3.2)
    thisuo2.add_element('O', 2)
    #thisuo2.add_element('Gd', 0.0007)
    thisuo2.set_density('g/cc', 10.341)
    uo2s.append(thisuo2)

materials = openmc.Materials(uo2s+waters)
"""






materials = openmc.Materials([uo2, water])
materials.export_to_xml()


L = pitch

fCylinders = [ openmc.ZCylinder(R=radius_fuel, x0=0.5*L+i*L, y0=0.5*L+j*L) for j in range(3) for i in range(3) ]

x1 = openmc.XPlane(x0=0.0, boundary_type='reflective')
x2 = openmc.XPlane(x0=1*L)
x3 = openmc.XPlane(x0=2*L)
x4 = openmc.XPlane(x0=3*L, boundary_type='reflective')

y1 = openmc.YPlane(y0=0.0, boundary_type='reflective')
y2 = openmc.YPlane(y0=1*L)
y3 = openmc.YPlane(y0=2*L)
y4 = openmc.YPlane(y0=3*L, boundary_type='reflective')

xP = [x1,x2,x3,x4]
yP = [y1,y2,y3,y4]



waterReg = +x1 & -x4 & +y1 & -y4 &                              \
           +fCylinders[0] &  +fCylinders[1] &  +fCylinders[2] & \
           +fCylinders[3] &  +fCylinders[4] &  +fCylinders[5] & \
           +fCylinders[6] &  +fCylinders[7] &  +fCylinders[8] 

fCells = [openmc.Cell(name='fuel'+str(i), fill=uo2, region=-fCylinders[i]) for i in range(9)]
#fCells = [openmc.Cell(name='fuel'+str(i), fill=uo2s[i], region=-fCylinders[i]) for i in range(9)]

mCells = []
count = 0
for j in range(3):
    for i in range(3):
        mRegion = waterReg & +xP[i] & -xP[i+1] & +yP[j] & -yP[j+1]
        mCells.append(openmc.Cell(70+count,'mod'+str(count),fill=water,region=mRegion))
        #mCells.append(openmc.Cell(70+count,'mod'+str(count),fill=waters[count],region=mRegion))
        count += 1



root = openmc.Universe(cells=(                                              \
    mCells[0], mCells[1], mCells[2], mCells[3], mCells[4], mCells[5], \
    mCells[6], mCells[7], mCells[8], fCells[0], fCells[1], fCells[2], \
    fCells[3], fCells[4], fCells[5], fCells[6], fCells[7], fCells[8] ))

geometry = openmc.Geometry(root)
geometry.export_to_xml()


# Settings

settings = openmc.Settings()
settings.batches = 100
settings.inactive = 25
settings.particles = 5000


space = openmc.stats.Box((0.0, 0.0, 0.0),(3.0*pitch, 3.0*pitch, 0))
settings.source = openmc.Source(space=space)
settings.export_to_xml()

groups = openmc.mgxs.EnergyGroups([0.0,0.625,20e6])

mgxs_lib = openmc.mgxs.Library(geometry)
mgxs_lib.energy_groups = groups

mgxs_lib.correction = 'P0'

mgxs_lib.mgxs_types = ('total','nu-transport','absorption','nu-fission', 'fission','consistent nu-scatter matrix','chi')

mgxs_lib.domain_type = 'cell'
mgxs_lib.domains = geometry.get_all_material_cells().values()
mgxs_lib.build_library()



tallies = openmc.Tallies()
mgxs_lib.add_to_tallies_file(tallies)
tallies.export_to_xml()


openmc.run()


sp = openmc.StatePoint('statepoint.100.h5')
mgxs_lib.load_from_statepoint(sp)

for i in range(9): mgxs_lib.domains[i].name   = 'fuel' + str(i)
for i in range(9): mgxs_lib.domains[9+i].name = 'mod'  + str(i)

mgxs_file = mgxs_lib.create_mg_library(xs_type='macro',xsdata_names=['fuel0','fuel1','fuel2','fuel3','fuel4','fuel5','fuel6','fuel7','fuel8','mod0','mod1','mod2','mod3','mod4','mod5','mod6','mod7','mod8'])


mgxs_file.export_to_hdf5()

fDatas = [ mgxs_file.get_by_name('fuel'+str(i)) for i in range(9) ]
mDatas = [ mgxs_file.get_by_name('mod' +str(i)) for i in range(9) ]

"""
print(fData.total)
print(fData.nu_fission)
print(fData.absorption)
print(fData.fission)
print(fData.chi)
print(fData.scatter_matrix)
print(fData.nu_fission)
"""

for fData in fDatas:
    assert(abs(fData.total[0][0]-                \
              (fData.absorption[0][0]+           \
               fData.scatter_matrix[0][0][0][0]+ \
               fData.scatter_matrix[0][0][1][0])) < 1e-12 )
for mData in mDatas:
    assert(abs(mData.total[0][0]-                \
              (mData.absorption[0][0]+           \
               mData.scatter_matrix[0][0][0][0]+ \
               mData.scatter_matrix[0][0][1][0])) < 1e-12 )




print(mDatas[0].scatter_matrix)

##################################################################
# PLOT
##################################################################
import sys
if (len(sys.argv) > 1):
    if (sys.argv[1] == 'plot'):

        p = openmc.Plot()
        p.filename = 'pinplot'
        p.width = (3*pitch, 3*pitch)
        p.pixels = (200, 200)
        p.color_by = 'material'
        p.colors = {uo2: 'yellow', water: 'blue'}
        p.origin = (3*pitch/2,3*pitch/2,0.0)


        plots = openmc.Plots([p])
        plots.export_to_xml()
        #get_ipython().system('cat plots.xml')
        openmc.plot_geometry()
        #get_ipython().system('convert pinplot.ppm pinplot.png')
        #from IPython.display import Image
        #Image("pinplot.png")

