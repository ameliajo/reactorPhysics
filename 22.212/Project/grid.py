

#get_ipython().magic('matplotlib inline')
import openmc
import openmc.mgxs as mgxs
import numpy as np


uo2 = openmc.Material(1,"uo2")
uo2.add_element('U', 1.0, enrichment=3.0)
uo2.add_element('O', 2.0)
uo2.set_density('g/cc', 10.0)

water = openmc.Material(3, "h2o")
water.add_element('H', 2.0)
water.add_element('O', 1.0)
water.set_density('g/cm3', 1.0)

mats = openmc.Materials([uo2, water])

mats.export_to_xml()


##################################################################
# GEOMETRY
##################################################################
L = 1.26
pitch = L*3

fuel_or_vec = [ openmc.ZCylinder(R=0.39, x0=(1+2*i)*L/2, y0=(1+2*j)*L/2) \
                for i in range(3) for j in range(3) ]


"""
fuel_vec = []
for i in range(9):
    fuel = openmc.Cell(i,'fuel'+str(i))
    fuel.fill = uo2
    fuel.region = -fuel_or_vec[i]
    fuel_vec.append(fuel)
"""
    
fuel_vec = [ openmc.Cell(i,'fuel'+str(i),fill=uo2,region=-fuel_or_vec[i]) for i in range(9)]


left   = openmc.XPlane(x0=0.0,   boundary_type='reflective')
right  = openmc.XPlane(x0=pitch, boundary_type='reflective')
bottom = openmc.YPlane(y0=0.0,   boundary_type='reflective')
top    = openmc.YPlane(y0=pitch, boundary_type='reflective')




water_region = +left & -right & +bottom & -top &                       \
               +fuel_or_vec[0] &  +fuel_or_vec[1] &  +fuel_or_vec[2] & \
               +fuel_or_vec[3] &  +fuel_or_vec[4] &  +fuel_or_vec[5] & \
               +fuel_or_vec[6] &  +fuel_or_vec[7] &  +fuel_or_vec[8] 



moderator = openmc.Cell(50, 'moderator')
moderator.fill = water
moderator.region = water_region


root = openmc.Universe(cells=(moderator,fuel_vec[0], fuel_vec[1], fuel_vec[2],\
                             fuel_vec[3], fuel_vec[4], fuel_vec[5],           \
                             fuel_vec[6], fuel_vec[7], fuel_vec[8] ))


geom = openmc.Geometry(root)
geom.export_to_xml()


##################################################################
# SETTINGS
##################################################################

settings = openmc.Settings()
settings.batches = 100
settings.inactive = 10
settings.particles = 10000
settings.output = {'tallies': True}

bounds = [0.0, 0.0, 0.0, pitch, pitch, pitch]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:])
settings.source = openmc.source.Source(space=uniform_dist)

settings.export_to_xml()



# Instantiate a 2-group EnergyGroups object
groups = mgxs.EnergyGroups()
groups.group_edges = np.array([0., 0.625, 20.0e6])



total_F = mgxs.TotalXS(domain=fuel_vec[0], groups=groups)
absorption_F = mgxs.AbsorptionXS(domain=fuel_vec[0], groups=groups)
scattering_F = mgxs.ScatterXS(domain=fuel_vec[0], groups=groups)

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
print(df)
df



##################################################################
# PLOT
##################################################################
"""
p = openmc.Plot()
p.filename = 'pinplot'
p.width = (pitch, pitch)
p.pixels = (200, 200)
p.color_by = 'material'
p.colors = {uo2: 'yellow', water: 'blue'}
p.origin = (pitch/2,pitch/2,0.0)




plots = openmc.Plots([p])
plots.export_to_xml()
#get_ipython().system('cat plots.xml')
openmc.plot_geometry()
#get_ipython().system('convert pinplot.ppm pinplot.png')
#from IPython.display import Image
#Image("pinplot.png")

"""
