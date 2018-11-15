
import openmc
import openmc.mgxs as mgxs
import numpy as np


uo2 = openmc.Material(1,"uo2")
uo2.add_element('U', 1.0, enrichment=3.0)
uo2.add_element('O', 2.0)
uo2.set_density('g/cc', 10.0)


zirconium = openmc.Material(2, "zirconium")
zirconium.add_element('Zr', 1.0)
zirconium.set_density('g/cm3', 6.6)

water = openmc.Material(3, "h2o")
water.add_nuclide('H1', 2.0)
water.add_nuclide('O16', 1.0)
water.set_density('g/cm3', 1.0)

mats = openmc.Materials([uo2, zirconium, water])


mats.export_to_xml()

##################################################################
# GEOMETRY
##################################################################
sideLen = 1.26
pitch = sideLen*3

fuel_or = openmc.ZCylinder(R=0.39,x0=sideLen/2,y0=sideLen/2)
clad_ir = openmc.ZCylinder(R=0.40,x0=sideLen/2,y0=sideLen/2)
clad_or = openmc.ZCylinder(R=0.46,x0=sideLen/2,y0=sideLen/2)

fuel_or2 = openmc.ZCylinder(R=0.39,x0=3*sideLen/2,y0=1*sideLen/2)
clad_ir2 = openmc.ZCylinder(R=0.40,x0=3*sideLen/2,y0=1*sideLen/2)
clad_or2 = openmc.ZCylinder(R=0.46,x0=3*sideLen/2,y0=1*sideLen/2)

fuel_or3 = openmc.ZCylinder(R=0.39,x0=5*sideLen/2,y0=1*sideLen/2)
clad_ir3 = openmc.ZCylinder(R=0.40,x0=5*sideLen/2,y0=1*sideLen/2)
clad_or3 = openmc.ZCylinder(R=0.46,x0=5*sideLen/2,y0=1*sideLen/2)

fuel_or4 = openmc.ZCylinder(R=0.39,x0=sideLen/2,y0=3*sideLen/2)
clad_ir4 = openmc.ZCylinder(R=0.40,x0=sideLen/2,y0=3*sideLen/2)
clad_or4 = openmc.ZCylinder(R=0.46,x0=sideLen/2,y0=3*sideLen/2)

fuel_or5 = openmc.ZCylinder(R=0.39,x0=3*sideLen/2,y0=3*sideLen/2)
clad_ir5 = openmc.ZCylinder(R=0.40,x0=3*sideLen/2,y0=3*sideLen/2)
clad_or5 = openmc.ZCylinder(R=0.46,x0=3*sideLen/2,y0=3*sideLen/2)

fuel_or6 = openmc.ZCylinder(R=0.39,x0=5*sideLen/2,y0=3*sideLen/2)
clad_ir6 = openmc.ZCylinder(R=0.40,x0=5*sideLen/2,y0=3*sideLen/2)
clad_or6 = openmc.ZCylinder(R=0.46,x0=5*sideLen/2,y0=3*sideLen/2)

fuel_or7 = openmc.ZCylinder(R=0.39,x0=sideLen/2,y0=5*sideLen/2)
clad_ir7 = openmc.ZCylinder(R=0.40,x0=sideLen/2,y0=5*sideLen/2)
clad_or7 = openmc.ZCylinder(R=0.46,x0=sideLen/2,y0=5*sideLen/2)

fuel_or8 = openmc.ZCylinder(R=0.39,x0=3*sideLen/2,y0=5*sideLen/2)
clad_ir8 = openmc.ZCylinder(R=0.40,x0=3*sideLen/2,y0=5*sideLen/2)
clad_or8 = openmc.ZCylinder(R=0.46,x0=3*sideLen/2,y0=5*sideLen/2)

fuel_or9 = openmc.ZCylinder(R=0.39,x0=5*sideLen/2,y0=5*sideLen/2)
clad_ir9 = openmc.ZCylinder(R=0.40,x0=5*sideLen/2,y0=5*sideLen/2)
clad_or9 = openmc.ZCylinder(R=0.46,x0=5*sideLen/2,y0=5*sideLen/2)



fuel_region = -fuel_or
gap_region = +fuel_or & -clad_ir
clad_region = +clad_ir & -clad_or

fuel_region2 = -fuel_or2
gap_region2 = +fuel_or2 & -clad_ir2
clad_region2 = +clad_ir2 & -clad_or2

fuel_region3 = -fuel_or3
gap_region3 = +fuel_or3 & -clad_ir3
clad_region3 = +clad_ir3 & -clad_or3

fuel_region4 = -fuel_or4
gap_region4 = +fuel_or4 & -clad_ir4
clad_region4 = +clad_ir4 & -clad_or4

fuel_region5 = -fuel_or5
gap_region5 = +fuel_or5 & -clad_ir5
clad_region5 = +clad_ir5 & -clad_or5

fuel_region6 = -fuel_or6
gap_region6 = +fuel_or6 & -clad_ir6
clad_region6 = +clad_ir6 & -clad_or6


fuel_region7 = -fuel_or7
gap_region7 = +fuel_or7 & -clad_ir7
clad_region7 = +clad_ir7 & -clad_or7

fuel_region8 = -fuel_or8
gap_region8 = +fuel_or8 & -clad_ir8
clad_region8 = +clad_ir8 & -clad_or8

fuel_region9 = -fuel_or9
gap_region9 = +fuel_or9 & -clad_ir9
clad_region9 = +clad_ir9 & -clad_or9



fuel = openmc.Cell(1, 'fuel')
fuel.fill = uo2
fuel.region = fuel_region

gap = openmc.Cell(2, 'air gap')
gap.region = gap_region

clad = openmc.Cell(3, 'clad')
clad.fill = zirconium
clad.region = clad_region


fuel2 = openmc.Cell(4, 'fuel2')
fuel2.fill = uo2
fuel2.region = fuel_region2

gap2 = openmc.Cell(5, 'air gap2')
gap2.region = gap_region2

clad2 = openmc.Cell(6, 'clad2')
clad2.fill = zirconium
clad2.region = clad_region2


fuel3 = openmc.Cell(7, 'fuel3')
fuel3.fill = uo2
fuel3.region = fuel_region3

gap3 = openmc.Cell(8, 'air gap3')
gap3.region = gap_region3

clad3 = openmc.Cell(9, 'clad3')
clad3.fill = zirconium
clad3.region = clad_region3



fuel4 = openmc.Cell(10, 'fuel4')
fuel4.fill = uo2
fuel4.region = fuel_region4

gap4 = openmc.Cell(11, 'air gap4')
gap4.region = gap_region4

clad4 = openmc.Cell(12, 'clad4')
clad4.fill = zirconium
clad4.region = clad_region4


fuel5 = openmc.Cell(13, 'fuel5')
fuel5.fill = uo2
fuel5.region = fuel_region5

gap5 = openmc.Cell(14, 'air gap5')
gap5.region = gap_region5

clad5 = openmc.Cell(15, 'clad5')
clad5.fill = zirconium
clad5.region = clad_region5


fuel6 = openmc.Cell(16, 'fuel6')
fuel6.fill = uo2
fuel6.region = fuel_region6

gap6 = openmc.Cell(17, 'air gap6')
gap6.region = gap_region6

clad6 = openmc.Cell(18, 'clad6')
clad6.fill = zirconium
clad6.region = clad_region6


fuel7 = openmc.Cell(19, 'fuel7')
fuel7.fill = uo2
fuel7.region = fuel_region7

gap7 = openmc.Cell(20, 'air gap7')
gap7.region = gap_region7

clad7 = openmc.Cell(21, 'clad7')
clad7.fill = zirconium
clad7.region = clad_region7


fuel8 = openmc.Cell(22, 'fuel8')
fuel8.fill = uo2
fuel8.region = fuel_region8

gap8 = openmc.Cell(23, 'air gap8')
gap8.region = gap_region8

clad8 = openmc.Cell(24, 'clad8')
clad8.fill = zirconium
clad8.region = clad_region8


fuel9 = openmc.Cell(25, 'fuel9')
fuel9.fill = uo2
fuel9.region = fuel_region9

gap9 = openmc.Cell(26, 'air gap9')
gap9.region = gap_region9

clad9 = openmc.Cell(27, 'clad9')
clad9.fill = zirconium
clad9.region = clad_region9



# left = openmc.XPlane(x0=-pitch/2, boundary_type='reflective')
# right = openmc.XPlane(x0=pitch/2, boundary_type='reflective')
# bottom = openmc.YPlane(y0=-pitch/2, boundary_type='reflective')
# top = openmc.YPlane(y0=pitch/2, boundary_type='reflective')

left = openmc.XPlane(x0=0.0, boundary_type='reflective')
right = openmc.XPlane(x0=pitch, boundary_type='reflective')
bottom = openmc.YPlane(y0=0.0, boundary_type='reflective')
top = openmc.YPlane(y0=pitch, boundary_type='reflective')


# left = openmc.XPlane(x0=-pitch/2, boundary_type='reflective')
# right = openmc.XPlane(x0=pitch/2, boundary_type='reflective')
# bottom = openmc.YPlane(y0=-pitch/2, boundary_type='reflective')
# top = openmc.YPlane(y0=pitch/2, boundary_type='reflective')



water_region = +left & -right & +bottom & -top &   \
               +clad_or & +clad_or2 & +clad_or3 &  \
               +clad_or4 & +clad_or5 & +clad_or6 & \
               +clad_or7 & +clad_or8 & +clad_or9

moderator = openmc.Cell(50, 'moderator')
moderator.fill = water
moderator.region = water_region

root = openmc.Universe(cells=(fuel, gap, clad, moderator,\
                              fuel2,gap2,clad2,\
                              fuel3,gap3,clad3,\
                              fuel4,gap4,clad4,\
                              fuel5,gap5,clad5,\
                              fuel6,gap6,clad6,\
                              fuel7,gap7,clad7,\
                              fuel8,gap8,clad8,\
                              fuel9,gap9,clad9))


geom = openmc.Geometry(root)
geom.export_to_xml()

##################################################################
# SETTINGS
##################################################################

# point = openmc.stats.Point((0, 0, 0))
# src = openmc.Source(space=point)



# Export to "settings.xml"
# settings_file.export_to_xml()

settings = openmc.Settings()
# settings.source = src
settings.batches = 100
settings.inactive = 10
settings.particles = 1000
settings.output = {'tallies': True}

bounds = [0.0, 0.0, 0.0, pitch, pitch, pitch]
# bounds = [0.1,0.1,0.1,0.2,0.2,0.2]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:])
settings.source = openmc.source.Source(space=uniform_dist)

settings.export_to_xml()

# Instantiate a 2-group EnergyGroups object
groups = mgxs.EnergyGroups()
groups.group_edges = np.array([0., 0.625, 20.0e6])

total_F = mgxs.TotalXS(domain=fuel, groups=groups)
absorption_F = mgxs.AbsorptionXS(domain=fuel, groups=groups)
scattering_F = mgxs.ScatterXS(domain=fuel, groups=groups)

total_M = mgxs.TotalXS(domain=moderator, groups=groups)
absorption_M = mgxs.AbsorptionXS(domain=moderator, groups=groups)
scattering_M = mgxs.ScatterXS(domain=moderator, groups=groups)

tallies_file = openmc.Tallies()

tallies_file += total_F.tallies.values()
tallies_file += absorption_F.tallies.values()
tallies_file += scattering_F.tallies.values()

tallies_file += total_M.tallies.values()
tallies_file += absorption_M.tallies.values()
tallies_file += scattering_M.tallies.values()

tallies_file.export_to_xml()

openmc.run()


sp = openmc.StatePoint('statepoint.100.h5')


total_F.load_from_statepoint(sp)
absorption_F.load_from_statepoint(sp)
scattering_F.load_from_statepoint(sp)

total_M.load_from_statepoint(sp)
absorption_M.load_from_statepoint(sp)
scattering_M.load_from_statepoint(sp)

df = total_F.get_pandas_dataframe()
df

##################################################################
# SETTINGS
##################################################################
p = openmc.Plot()
p.filename = 'pinplot'
p.width = (pitch, pitch)
p.pixels = (200, 200)
p.color_by = 'material'
p.colors = {uo2: 'yellow', water: 'blue'}
p.origin = (pitch/2,pitch/2,0.0)


plots = openmc.Plots([p])
plots.export_to_xml()
#!cat plots.xml
openmc.plot_geometry()
#!convert pinplot.ppm pinplot.png
#from IPython.display import Image
#Image("pinplot.png")




