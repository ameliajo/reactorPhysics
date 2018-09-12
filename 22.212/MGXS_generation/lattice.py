import openmc
import matplotlib


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
    fuel.set_density('g/cc',m1.density)
    return fuel

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

pitch = 1.26

m1 = material(1); m2 = material(2); m3 = material(3);
m4 = material(4); m5 = material(5); m6 = material(6); 
m7 = material(7); m8 = material(8); m9 = material(9);

uo2_1 = makeFuel(m1); uo2_2 = makeFuel(m2); uo2_3 = makeFuel(m3); 
uo2_4 = makeFuel(m4); uo2_5 = makeFuel(m5); uo2_6 = makeFuel(m6); 
uo2_7 = makeFuel(m7); uo2_8 = makeFuel(m8); uo2_9 = makeFuel(m9);

water = openmc.Material(10,"water")
water.add_element('O',1.0)
water.add_element('H',2.0)
water.set_density('g/cc',1.0)
water.add_s_alpha_beta('c_H_in_H2O')

mats = openmc.Materials([uo2_1,uo2_2,uo2_3,uo2_4,uo2_5,uo2_6,uo2_7,uo2_8,uo2_9,water])
mats.export_to_xml()

#Define boundaries
fuel_or = openmc.ZCylinder(R=0.39218)

left  = openmc.XPlane(x0=0.0,name="left")
right = openmc.XPlane(x0=3.0*pitch,name="right")
down  = openmc.YPlane(y0=0.0,name="down")
up    = openmc.YPlane(y0=3.0*pitch,name="up")

left.boundary_type  = 'reflective'; right.boundary_type = 'reflective';
down.boundary_type  = 'reflective'; up.boundary_type    = 'reflective';

outer_box_region = +left & -right & +down & -up
outer_box_cell = openmc.Cell(name="OuterBox")
outer_box_cell.region = outer_box_region 
 

fuel_region1  = -fuel_or; fuel_region2  = -fuel_or; fuel_region3  = -fuel_or; 
fuel_region4  = -fuel_or; fuel_region5  = -fuel_or; fuel_region6  = -fuel_or; 
fuel_region7  = -fuel_or; fuel_region8  = -fuel_or; fuel_region9  = -fuel_or;

water_region1 = +fuel_or; water_region2 = +fuel_or; water_region3 = +fuel_or; 
water_region4 = +fuel_or; water_region5 = +fuel_or; water_region6 = +fuel_or; 
water_region7 = +fuel_or; water_region8 = +fuel_or; water_region9 = +fuel_or;

fuel_1 = makeFuelCell(uo2_1,fuel_region1,m1)
fuel_2 = makeFuelCell(uo2_2,fuel_region2,m2)
fuel_3 = makeFuelCell(uo2_3,fuel_region3,m3)
fuel_4 = makeFuelCell(uo2_4,fuel_region4,m4)
fuel_5 = makeFuelCell(uo2_5,fuel_region5,m5)
fuel_6 = makeFuelCell(uo2_6,fuel_region6,m6)
fuel_7 = makeFuelCell(uo2_7,fuel_region7,m7)
fuel_8 = makeFuelCell(uo2_8,fuel_region8,m8)
fuel_9 = makeFuelCell(uo2_9,fuel_region9,m9)

moderator_1 = makeWaterCell(water,water_region1,m1)
moderator_2 = makeWaterCell(water,water_region2,m2)
moderator_3 = makeWaterCell(water,water_region3,m3)
moderator_4 = makeWaterCell(water,water_region4,m4)
moderator_5 = makeWaterCell(water,water_region5,m5)
moderator_6 = makeWaterCell(water,water_region6,m6)
moderator_7 = makeWaterCell(water,water_region7,m7)
moderator_8 = makeWaterCell(water,water_region8,m8)
moderator_9 = makeWaterCell(water,water_region9,m9)

u_1 = openmc.Universe(cells=[fuel_1, moderator_1])
u_2 = openmc.Universe(cells=[fuel_2, moderator_2])
u_3 = openmc.Universe(cells=[fuel_3, moderator_3])
u_4 = openmc.Universe(cells=[fuel_4, moderator_4])
u_5 = openmc.Universe(cells=[fuel_5, moderator_5])
u_6 = openmc.Universe(cells=[fuel_6, moderator_6])
u_7 = openmc.Universe(cells=[fuel_7, moderator_7])
u_8 = openmc.Universe(cells=[fuel_8, moderator_8])
u_9 = openmc.Universe(cells=[fuel_9, moderator_9])


lattice = openmc.RectLattice()
lattice.lower_left = (0.0, 0.0)
lattice.pitch = (pitch, pitch)
lattice.universes = [[u_1, u_2, u_3],[u_4, u_5, u_6],[u_7, u_8, u_9]]

settings = openmc.Settings()
full = 3*pitch
bounds = [0,0,0,full,full,full]
uniform_dist = openmc.stats.Box(bounds[:3],bounds[3:],only_fissionable=True)
settings.source = openmc.source.Source(space=uniform_dist)
settings.batches = 50
settings.inactive = 10
settings.particles = 1000
settings.export_to_xml()

root = openmc.Universe()
root.add_cell(outer_box_cell)
outer_box_cell.fill = lattice
geom = openmc.Geometry(root)
geom.export_to_xml()

plot = openmc.Plot(plot_id=1)
plot.origin = [3, 3, 0]
plot.width = [8, 8]
plot.pixels = [1000, 1000]
plot.color_by = 'material'

# Instantiate a Plots collection and export to XML
plot_file = openmc.Plots([plot])
plot_file.export_to_xml()
plot = openmc.plot_geometry()

cell_filter = openmc.CellFilter([fuel_1, moderator_1])

energy_filter = openmc.EnergyFilter([0., 4.0, 20.0e6])
t = openmc.Tally(1)
t.filters = [cell_filter, energy_filter]
# these are the main reaction rates you should need
t.scores = ['absorption','nu-fission','fission']
tallies = openmc.Tallies([t])
tallies.export_to_xml()
openmc.run()


