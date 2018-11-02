import numpy as np
import matplotlib.pyplot as plt
plt.style.use('seaborn-dark')

import openmc
import openmc.mgxs as mgxs
import openmc.data

#%matplotlib inline
fuel = openmc.Material(name='Fuel')
fuel.set_density('g/cm3', 10.31341)
fuel.add_nuclide('U235', 0.03)
fuel.add_nuclide('U238', 0.97)
fuel.add_nuclide('O16', 2.0)

water = openmc.Material(name='Water')
water.set_density('g/cm3', 0.740582)
water.add_nuclide('H1', 2.0)
water.add_nuclide('O16', 1.0)
water.add_s_alpha_beta('c_H_in_H2O')

# Instantiate a Materials collection
materials_file = openmc.Materials([fuel, water])

# Export to "materials.xml"
materials_file.export_to_xml()

# Create cylinders for the fuel and clad
fuel_outer_radius = openmc.ZCylinder(x0=0.0, y0=0.0, R=0.39218)

# Create a Universe to encapsulate a fuel pin
pin_cell_universe = openmc.Universe(name='Fuel Pin')

# Create fuel Cell
fuel_cell = openmc.Cell(name='Fuel')
fuel_cell.fill = fuel
fuel_cell.region = -fuel_outer_radius
pin_cell_universe.add_cell(fuel_cell)

# Create a clad Cell
clad_cell = openmc.Cell(name='Clad')
# clad_cell.fill = fuel
# clad_cell.region = +fuel_outer_radius & -clad_outer_radius
# pin_cell_universe.add_cell(clad_cell)

# Create a moderator Cell
moderator_cell = openmc.Cell(name='Moderator')
moderator_cell.fill = water
moderator_cell.region = +fuel_outer_radius
pin_cell_universe.add_cell(moderator_cell)

# Create root Cell
root_cell = openmc.Cell(name='root cell')
box = openmc.get_rectangular_prism(width=1.26, height=1.26,
                                   boundary_type='reflective')
root_cell.region = box
root_cell.fill = pin_cell_universe

# Create root Universe
root_universe = openmc.Universe(universe_id=0, name='root universe')
root_universe.add_cell(root_cell)

# Create Geometry and set root Universe
openmc_geometry = openmc.Geometry(root_universe)

# Export to "geometry.xml"
openmc_geometry.export_to_xml()

# Instantiate a Settings object
settings_file = openmc.Settings()
settings_file.batches = 50                # 50
settings_file.inactive = 10               # 10
settings_file.particles = 1000             # 10000
settings_file.output = {'tallies': True}

# Create an initial uniform spatial source distribution over fissionable zones
bounds = [-0.63, -0.63, -0.63, 0.63, 0.63, 0.63]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
settings_file.source = openmc.source.Source(space=uniform_dist)

# Activate tally precision triggers
settings_file.trigger_active = True
settings_file.trigger_max_batches = settings_file.batches * 4

# Export to "settings.xml"
settings_file.export_to_xml()


# Instantiate a "coarse" 2-group EnergyGroups object
coarse_groups = mgxs.EnergyGroups([0., 0.625, 20.0e6])

# Instantiate a "fine" 8-group EnergyGroups object
fine_groups = mgxs.EnergyGroups([0., 0.058, 0.14, 0.28,
                                 0.625, 4.0, 5.53e3, 821.0e3, 20.0e6])


# Extract all Cells filled by Materials
openmc_cells = openmc_geometry.get_all_material_cells().values()

# Create dictionary to store multi-group cross sections for all cells
xs_library = {}

# Instantiate 8-group cross sections for each cell
for cell in openmc_cells:
    xs_library[cell.id] = {}
    xs_library[cell.id]['transport']  = mgxs.TransportXS(groups=fine_groups)
    xs_library[cell.id]['fission']    = mgxs.FissionXS(groups=fine_groups)
    xs_library[cell.id]['nu-fission'] = mgxs.FissionXS(groups=fine_groups, nu=True)
    xs_library[cell.id]['nu-scatter'] = mgxs.ScatterMatrixXS(groups=fine_groups, nu=True)
    xs_library[cell.id]['chi']        = mgxs.Chi(groups=fine_groups)



# Create a tally trigger for +/- 0.01 on each tally used to compute the multi-group cross sections
tally_trigger = openmc.Trigger('std_dev', 1E-2)

# Add the tally trigger to each of the multi-group cross section tallies
for cell in openmc_cells:
    for mgxs_type in xs_library[cell.id]:
        xs_library[cell.id][mgxs_type].tally_trigger = tally_trigger



# Instantiate an empty Tallies object
tallies_file = openmc.Tallies()

# Iterate over all cells and cross section types
for cell in openmc_cells:
    for rxn_type in xs_library[cell.id]:

        # Set the cross sections domain to the cell
        xs_library[cell.id][rxn_type].domain = cell
        
        # Tally cross sections by nuclide
        xs_library[cell.id][rxn_type].by_nuclide = True
                
        # Add OpenMC tallies to the tallies file for XML generation
        for tally in xs_library[cell.id][rxn_type].tallies.values():
            tallies_file.append(tally, merge=True)

# Export to "tallies.xml"
tallies_file.export_to_xml()

# Run OpenMC
openmc.run()

# Load the last statepoint file
sp = openmc.StatePoint('statepoint.200.h5')

# Iterate over all cells and cross section types
for cell in openmc_cells:
    for rxn_type in xs_library[cell.id]:
        xs_library[cell.id][rxn_type].load_from_statepoint(sp)

# nufission = xs_library[fuel_cell.id]['nu-fission']
# nufission.print_xs(xs_type='micro', nuclides=['U235', 'U238'])


# nufission = xs_library[fuel_cell.id]['nu-fission']
# nufission.print_xs(xs_type='macro', nuclides='sum')

nufission = xs_library[moderator_cell.id]['nu-fission']
df = nufission.get_pandas_dataframe(xs_type='macro')
print("Nu-Fission Macro in Moderator")
print(df)
print()

nufission = xs_library[fuel_cell.id]['nu-fission']
df = nufission.get_pandas_dataframe(xs_type='macro')
print("Nu-Fission Macro in Fuel")
print(df)
print()





"""
nuscatter = xs_library[moderator_cell.id]['nu-scatter']
df = nuscatter.get_pandas_dataframe(xs_type='micro')
print("Nu-Scatter")
print(df.head(10))
print()

nuscatter = xs_library[moderator_cell.id]['nu-scatter']
df = nuscatter.get_pandas_dataframe(xs_type='micro')
print(df.head(10))
"""

#df = condensed_xs.get_pandas_dataframe(xs_type='micro')
#print(df)

# Create a figure of the U-235 continuous-energy fission cross section
fig = openmc.plot_xs('U235', ['fission'])

# Get the axis to use for plotting the MGXS
ax = fig.gca()

# Extract energy group bounds and MGXS values to plot
fission = xs_library[fuel_cell.id]['fission']
energy_groups = fission.energy_groups
x = energy_groups.group_edges
y = fission.get_xs(nuclides=['U235'], order_groups='decreasing', xs_type='micro')
y = np.squeeze(y)

# Fix low energy bound
x[0] = 1.e-5

# Extend the mgxs values array for matplotlib's step plot
y = np.insert(y, 0, y[0])

# Create a step plot for the MGXS
ax.plot(x, y, drawstyle='steps', color='r', linewidth=3)

ax.set_title('U-235 Fission Cross Section')
ax.legend(['Continuous', 'Multi-Group'])
ax.set_xlim((x.min(), x.max()))

plt.show()





