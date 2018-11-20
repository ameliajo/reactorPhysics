from math import pi
import openmc
import openmc.mgxs
import openmc.model
import numpy as np


radius_fuel = 0.3922
radius_gap = 0.4001
radius_clad = 0.4572
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


materials = openmc.Materials([uo2, water])
materials.export_to_xml()


L = pitch

x1 = openmc.XPlane(x0=0.0, boundary_type='reflective')
x2 = openmc.XPlane(x0=1*L)
x3 = openmc.XPlane(x0=2*L)
x4 = openmc.XPlane(x0=3*L, boundary_type='reflective')

y1 = openmc.YPlane(y0=0.0, boundary_type='reflective')
y2 = openmc.YPlane(y0=1*L)
y3 = openmc.YPlane(y0=2*L)
y4 = openmc.YPlane(y0=3*L, boundary_type='reflective')





rf0 = openmc.ZCylinder(R=radius_fuel,x0=pitch*0.5,y0=pitch*0.5)
rf1 = openmc.ZCylinder(R=radius_fuel,x0=pitch*0.5,y0=pitch*1.5)
rf2 = openmc.ZCylinder(R=radius_fuel,x0=pitch*0.5,y0=pitch*2.5)
rf3 = openmc.ZCylinder(R=radius_fuel,x0=pitch*1.5,y0=pitch*0.5)
rf4 = openmc.ZCylinder(R=radius_fuel,x0=pitch*1.5,y0=pitch*1.5)
rf5 = openmc.ZCylinder(R=radius_fuel,x0=pitch*1.5,y0=pitch*2.5)
rf6 = openmc.ZCylinder(R=radius_fuel,x0=pitch*2.5,y0=pitch*0.5)
rf7 = openmc.ZCylinder(R=radius_fuel,x0=pitch*2.5,y0=pitch*1.5)
rf8 = openmc.ZCylinder(R=radius_fuel,x0=pitch*2.5,y0=pitch*2.5)

waterReg = +rf0 & +rf1 & +rf2 & +rf3 & +rf4 & +rf5 & +rf6 & +rf7 & +rf8 & +x1 & -x4 & +y1 & -y4

fuel0 = openmc.Cell(name='fuel0', fill=uo2, region=-rf0)
fuel1 = openmc.Cell(name='fuel1', fill=uo2, region=-rf1)
fuel2 = openmc.Cell(name='fuel2', fill=uo2, region=-rf2)
fuel3 = openmc.Cell(name='fuel3', fill=uo2, region=-rf3)
fuel4 = openmc.Cell(name='fuel4', fill=uo2, region=-rf4)
fuel5 = openmc.Cell(name='fuel5', fill=uo2, region=-rf5)
fuel6 = openmc.Cell(name='fuel6', fill=uo2, region=-rf6)
fuel7 = openmc.Cell(name='fuel7', fill=uo2, region=-rf7)
fuel8 = openmc.Cell(name='fuel8', fill=uo2, region=-rf8)

mod = openmc.Cell(name='moderator', fill=water, region=waterReg)
root = openmc.Universe(cells=(fuel0,fuel1,fuel2,fuel3,fuel4,fuel5,fuel6,fuel7,fuel8,mod))
geometry = openmc.Geometry(root)
geometry.export_to_xml()


# Settings

settings = openmc.Settings()

settings.batches = 100
settings.inactive = 25
settings.particles = 5000


# Set the initial source to a flat distribution born only in the fuel.
space = openmc.stats.Box((0.0, 0.0, 0.0),
     (3*pitch, 3*pitch, 0))
settings.source = openmc.Source(space=space)
settings.export_to_xml()




# # MGXS Tallies
groups = openmc.mgxs.EnergyGroups([0.0,0.625,20e6])

# Instantiate an MGXS library.
mgxs_lib = openmc.mgxs.Library(geometry)
mgxs_lib.energy_groups = groups

# Don't apply any anisotropic scattering corrections.
mgxs_lib.correction = 'P0'

# Set the desired MGXS data.
mgxs_lib.mgxs_types = ('total','nu-transport','absorption','nu-fission', 'fission','consistent nu-scatter matrix','chi')

# Define the domain and build the library.
mgxs_lib.domain_type = 'cell'
mgxs_lib.domains = geometry.get_all_material_cells().values()
#mgxs_lib.domain_type = 'universe'
#mgxs_lib.domains = geometry.get_all_material_universes().values()
mgxs_lib.build_library()


# Add the tallies.
tallies = openmc.Tallies()
mgxs_lib.add_to_tallies_file(tallies)

tallies.export_to_xml()

openmc.run()


# Load the statepoint and the MGXS results.
sp = openmc.StatePoint('statepoint.100.h5')
mgxs_lib.load_from_statepoint(sp)


# Pick out the fuel and moderator domains.
#fuel0 = mgxs_lib.domains[0]
fuel0 = mgxs_lib.domains[0]
fuel0.name = 'fuel0'
fuel1 = mgxs_lib.domains[1]
fuel1.name = 'fuel1'
fuel2 = mgxs_lib.domains[2]
fuel2.name = 'fuel2'
fuel3 = mgxs_lib.domains[3]
fuel3.name = 'fuel3'
fuel4 = mgxs_lib.domains[4]
fuel4.name = 'fuel4'
fuel5 = mgxs_lib.domains[5]
fuel5.name = 'fuel5'
fuel6 = mgxs_lib.domains[6]
fuel6.name = 'fuel6'
fuel7 = mgxs_lib.domains[7]
fuel7.name = 'fuel7'
fuel8 = mgxs_lib.domains[8]
fuel8.name = 'fuel8'


moderator = mgxs_lib.domains[9]
moderator.name = 'moderator'
#assert moderator.name == 'moderator'


mgxs_file = mgxs_lib.create_mg_library(xs_type='macro',xsdata_names=['fuel0','fuel1','fuel2','fuel3','fuel4','fuel5','fuel6','fuel7','fuel8','moderator'])
#mgxs_file = mgxs_lib.create_mg_library(xs_type='macro',xsdata_names=['fuel0'])


mgxs_file.export_to_hdf5()


f0Data = mgxs_file.get_by_name('fuel0')
f1Data = mgxs_file.get_by_name('fuel1')
f2Data = mgxs_file.get_by_name('fuel2')
f3Data = mgxs_file.get_by_name('fuel3')
f4Data = mgxs_file.get_by_name('fuel4')
f5Data = mgxs_file.get_by_name('fuel5')
f6Data = mgxs_file.get_by_name('fuel6')
f7Data = mgxs_file.get_by_name('fuel7')
f8Data = mgxs_file.get_by_name('fuel8')
mData = mgxs_file.get_by_name('moderator')

"""
print(fData.total)
print(fData.nu_fission)
print(fData.absorption)
print(fData.fission)
print(fData.chi)
print(fData.scatter_matrix)
print(fData.nu_fission)
"""

fDatas = [f0Data,f1Data,f2Data,f3Data,f4Data,f5Data,f6Data,f7Data,f8Data]
for fData in fDatas:
    assert(abs(fData.total[0][0]-                \
              (fData.absorption[0][0]+           \
               fData.scatter_matrix[0][0][0][0]+ \
               fData.scatter_matrix[0][0][1][0])) < 1e-12 )



"""
print(f0Data.total[0][0])
print((f0Data.absorption[0][0]+           \
       f0Data.scatter_matrix[0][0][0][0]+ \
       f0Data.scatter_matrix[0][0][1][0]))

assert(abs(f0Data.total[0][0]-                \
          (f0Data.absorption[0][0]+           \
           f0Data.scatter_matrix[0][0][0][0]+ \
           f0Data.scatter_matrix[0][0][1][0])) < 1e-12 )


print(f1Data.total[0][0])
print((f1Data.absorption[0][0]+           \
       f1Data.scatter_matrix[0][0][0][0]+ \
       f1Data.scatter_matrix[0][0][1][0]))

assert(abs(f1Data.total[0][0]-                \
          (f1Data.absorption[0][0]+           \
           f1Data.scatter_matrix[0][0][0][0]+ \
           f1Data.scatter_matrix[0][0][1][0])) < 1e-12 )



print(mData.total[0][0])
print((mData.absorption[0][0]+           \
       mData.scatter_matrix[0][0][0][0]+ \
       mData.scatter_matrix[0][0][1][0]))

assert(abs(mData.total[0][0]-                \
          (mData.absorption[0][0]+           \
           mData.scatter_matrix[0][0][0][0]+ \
           mData.scatter_matrix[0][0][1][0])) < 1e-12 )

print(mData.scatter_matrix)
"""

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

