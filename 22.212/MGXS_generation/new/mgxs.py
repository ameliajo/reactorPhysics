import numpy as np
import matplotlib.pyplot as plt
plt.style.use('seaborn-dark')
import openmc
import openmc.mgxs as mgxs
import openmc.data
#%matplotlib inline



# ----------------------------------------------------------------------------
# Materials
# ----------------------------------------------------------------------------

fuel = openmc.Material(name='Fuel')
fuel.add_element('U', 1.0, enrichment=3.0)
fuel.add_element('O', 2.0)
fuel.set_density('g/cc', 10.31341)

water = openmc.Material(name='Water')
water.set_density('g/cm3', 0.740582)
water.add_element('H', 2.0)
water.add_element('O', 1.0)
water.add_s_alpha_beta('c_H_in_H2O')

materials_file = openmc.Materials([fuel, water])
materials_file.export_to_xml()



# ----------------------------------------------------------------------------
# Geometry 
# ----------------------------------------------------------------------------

fuel_outer_radius = openmc.ZCylinder(x0=0.0, y0=0.0, R=0.39218)
pin_cell_universe = openmc.Universe(name='Fuel Pin')

fuel_cell = openmc.Cell(name='Fuel')
fuel_cell.fill = fuel
fuel_cell.region = -fuel_outer_radius
pin_cell_universe.add_cell(fuel_cell)

moderator_cell = openmc.Cell(name='Moderator')
moderator_cell.fill = water
moderator_cell.region = +fuel_outer_radius
pin_cell_universe.add_cell(moderator_cell)

root_cell = openmc.Cell(name='root cell')
box = openmc.get_rectangular_prism(width=1.26, height=1.26,
                                   boundary_type='reflective')
root_cell.region = box
root_cell.fill = pin_cell_universe

root_universe = openmc.Universe(universe_id=0, name='root universe')
root_universe.add_cell(root_cell)

openmc_geometry = openmc.Geometry(root_universe)
openmc_geometry.export_to_xml()



# ----------------------------------------------------------------------------
# Settings 
# ----------------------------------------------------------------------------

settings_file = openmc.Settings()
settings_file.batches = 50                # 50
settings_file.inactive = 10               # 10
settings_file.particles = 1000             # 10000
settings_file.output = {'tallies': True}

bounds = [-0.63, -0.63, -0.63, 0.63, 0.63, 0.63]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
settings_file.source = openmc.source.Source(space=uniform_dist)

settings_file.trigger_active = True
settings_file.trigger_max_batches = settings_file.batches * 4
settings_file.export_to_xml()


coarse_groups = mgxs.EnergyGroups([0., 0.625, 20e6])
fine_groups = mgxs.EnergyGroups([0., 0.058, 0.14, 0.28, 0.625, 4, 5.53e3, 821e3, 20e6])


openmc_cells = openmc_geometry.get_all_material_cells().values()

xsLib = {}
for cell in openmc_cells:
    xsLib[cell.id] = {}
    xsLib[cell.id]['transport']  = mgxs.TransportXS(    groups=fine_groups)
    xsLib[cell.id]['fission']    = mgxs.FissionXS(      groups=fine_groups)
    xsLib[cell.id]['nu-fission'] = mgxs.FissionXS(      groups=fine_groups, nu=True)
    xsLib[cell.id]['nu-scatter'] = mgxs.ScatterMatrixXS(groups=fine_groups, nu=True)
    xsLib[cell.id]['chi']        = mgxs.Chi(            groups=fine_groups)



tally_trigger = openmc.Trigger('std_dev', 1E-2)

for cell in openmc_cells:
    for mgxs_type in xsLib[cell.id]:
        xsLib[cell.id][mgxs_type].tally_trigger = tally_trigger



tallies_file = openmc.Tallies()

# Iterate over all cells and cross section types
for cell in openmc_cells:
    for rxn_type in xsLib[cell.id]:

        xsLib[cell.id][rxn_type].domain = cell
        xsLib[cell.id][rxn_type].by_nuclide = True
                
        for tally in xsLib[cell.id][rxn_type].tallies.values():
            tallies_file.append(tally, merge=True)

tallies_file.export_to_xml()
openmc.run()
sp = openmc.StatePoint('statepoint.200.h5')

for cell in openmc_cells:
    for rxn_type in xsLib[cell.id]:
        xsLib[cell.id][rxn_type].load_from_statepoint(sp)

#nufission = xsLib[fuel_cell.id]['nu-fission']
#nufission.print_xs(xs_type='micro', nuclides=['U235', 'U238'])

nufission = xsLib[fuel_cell.id]['nu-fission']
df = nufission.get_pandas_dataframe(xs_type='macro')
print("Nu-Fission Macro in Fuel")
print(df)
print()






"""
fig = openmc.plot_xs('U235', ['fission'])
ax = fig.gca()

fission = xsLib[fuel_cell.id]['fission']
energy_groups = fission.energy_groups
x = energy_groups.group_edges
y = fission.get_xs(nuclides=['U235'], order_groups='decreasing', xs_type='micro')
y = np.squeeze(y)

x[0] = 1.e-5

y = np.insert(y, 0, y[0])

ax.plot(x, y, drawstyle='steps', color='r', linewidth=3)

ax.set_title('U-235 Fission Cross Section')
ax.legend(['Continuous', 'Multi-Group'])
ax.set_xlim((x.min(), x.max()))

plt.show()
"""





