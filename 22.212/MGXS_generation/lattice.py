import openmc
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import openmc.mgxs as mgxs


class material:
    def __init__(self,idVal,density=10.29769):
        self.id   = idVal
        self.density = density
    def __str__(self):
        return "\nMaterial #%s\nTemp         %s\nEnrichment   %s\nDensity      %s\n" % (self.id,self.temp,self.enrichment,self.density)

def makeFuel(m):
    fuel = openmc.Material(m.id,"fuel")
    fuel.add_element('U',1.0,enrichment=3.035)
    fuel.add_element('O',2.0)
    fuel.set_density('g/cc',m.density)
    return fuel

def makeWater():
    water = openmc.Material(10,"water")
    water.add_element('O',1.0)
    water.add_element('H',2.0)
    water.set_density('g/cc',1.0)
    water.add_s_alpha_beta('c_H_in_H2O')
    return water

def makeFuelCell(uo2,fuel_region,m):
    fuelCell = openmc.Cell(name='fuel')
    fuelCell.fill = uo2
    fuelCell.region = fuel_region
    return fuelCell

def makeWaterCell(water,water_region,m):
    moderator = openmc.Cell(name='moderator')
    moderator.fill = water
    moderator.region = water_region 
    return moderator

def plotGeom():
    plot = openmc.Plot(plot_id=1)
    plot.origin = [3, 3, 0]
    plot.width = [8, 8]
    plot.pixels = [1000, 1000]
    plot.color_by = 'material'
    plot_file = openmc.Plots([plot])
    plot_file.export_to_xml()
    plot = openmc.plot_geometry()



# ---------------------------------------------------------------------------
# Materials
# ---------------------------------------------------------------------------
m   = [0]*9
uo2 = [0]*9
for i in range(9):
    m[i]   = material(i+1)  # Class object,  easier to change individual pins
    uo2[i] = makeFuel(m[i]) # Define each fuel cell

water = makeWater() # Define water cell

mats = openmc.Materials([uo2[0],uo2[1],uo2[2],uo2[3],uo2[4],uo2[5],uo2[6],
                         uo2[7],uo2[8],water])
mats.export_to_xml()



# ---------------------------------------------------------------------------
# Settings 
# ---------------------------------------------------------------------------
pitch = 1.26
settings = openmc.Settings()
bounds = [0,0,0,3.0*pitch,3.0*pitch,3.0*pitch]
uniform_dist = openmc.stats.Box(bounds[:3],bounds[3:],only_fissionable=True)
settings.source = openmc.source.Source(space=uniform_dist)
#settings.batches   = 50; settings.inactive  = 10; settings.particles = 1000
settings.batches   = 100
settings.inactive  = 10
settings.particles = 10000
settings.export_to_xml()


# ---------------------------------------------------------------------------
# Geometry 
# ---------------------------------------------------------------------------

fuel_or = openmc.ZCylinder(R=0.39218)

L = openmc.XPlane(x0=0.0,      name="left" )
R = openmc.XPlane(x0=3.0*pitch,name="right")
D = openmc.YPlane(y0=0.0,      name="down" )
U = openmc.YPlane(y0=3.0*pitch,name="up"   )

L.boundary_type = 'reflective'; R.boundary_type = 'reflective';
D.boundary_type = 'reflective'; U.boundary_type = 'reflective';

outer_box_cell  = openmc.Cell(name="OuterBox")
outer_box_cell.region = +L & -R & +D & -U
 

fuel         = [0]*9
moderator    = [0]*9
univ         = [0]*9

for i in range(9):
    fuel[i]      = makeFuelCell(uo2[i],-fuel_or,m[i])
    moderator[i] = makeWaterCell(water,+fuel_or,m[i])
    univ[i]      = openmc.Universe(cells=[fuel[i], moderator[i]])


lattice = openmc.RectLattice()
lattice.lower_left = (0.0, 0.0)
lattice.pitch = (pitch, pitch)
lattice.universes = [[univ[0],univ[1],univ[2]],[univ[3],univ[4],univ[5]],[univ[6],univ[7],univ[8]]]

root = openmc.Universe()
root.add_cell(outer_box_cell)
outer_box_cell.fill = lattice
geom = openmc.Geometry(root)
geom.export_to_xml()


# ---------------------------------------------------------------------------
# Tallies
# ---------------------------------------------------------------------------


# Instantiate a 2-group EnergyGroups object
groups = mgxs.EnergyGroups()
groups.group_edges = np.array([0.0, 0.058, 0.14, 0.28, 0.625, 4.0, 10.0, 40.0, 5530.0, 821e3, 20e6])
# Instantiate a few different sections

total = mgxs.TotalXS(domain=outer_box_cell, groups=groups)
absorption = mgxs.AbsorptionXS(domain=outer_box_cell, groups=groups)
scattering = mgxs.ScatterXS(domain=outer_box_cell, groups=groups)


#plotGeom()
openmc.run()


# Load the last statepoint file
sp = openmc.StatePoint('statepoint.50.h5')
# Load the tallies from the statepoint into each MGXS object
total.load_from_statepoint(sp)
absorption.load_from_statepoint(sp)
scattering.load_from_statepoint(sp)
total.print_xs()
absorption.print_xs()
scattering.print_xs()














