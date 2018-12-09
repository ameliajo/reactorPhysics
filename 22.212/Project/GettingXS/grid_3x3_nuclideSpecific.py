from math import pi
import openmc
import openmc.mgxs
import openmc.model
import numpy as np

import sys
#sys.stdout = open('output','w')



radius_fuel = 0.39128
pitch = 1.26

# Basic materials
uo2 = openmc.Material(name='fuel')
uo2.add_element('U', 1, enrichment=3.2)
uo2.add_element('O', 2)
#uo2.add_element('Gd', 0.0007)
uo2.set_density('g/cc', 10.341)
x = uo2.get_nuclide_atom_densities()
#for d in x:
#    print(d,x[d],x[d][1]*1E24)




water = openmc.Material(3, "h2o")
water.add_element('H', 2.0)
water.add_element('O', 1.0)
water.set_density('g/cm3', 1.0)


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

mCells = []
count = 0
for j in range(3):
    for i in range(3):
        mReg = waterReg & +xP[i] & -xP[i+1] & +yP[j] & -yP[j+1]
        mCells.append(openmc.Cell(70+count,'mod'+str(count),fill=water,region=mReg))
        count += 1


root = openmc.Universe(cells=(                                        \
    mCells[0], mCells[1], mCells[2], mCells[3], mCells[4], mCells[5], \
    mCells[6], mCells[7], mCells[8], fCells[0], fCells[1], fCells[2], \
    fCells[3], fCells[4], fCells[5], fCells[6], fCells[7], fCells[8] ))

geometry = openmc.Geometry(root)
geometry.export_to_xml()


# Settings

settings = openmc.Settings()
settings.batches = 100
settings.inactive = 25
settings.particles = 500



space = openmc.stats.Box((0.0, 0.0, 0.0),(3.0*pitch, 3.0*pitch, 0))
settings.source = openmc.Source(space=space)
settings.export_to_xml()

groups = [0.0, 4, 10, 20e6]
groups = [0.0, 0.058, 4, 10, 20e6]
groups = [0.0, 20e6]
groups = [0.0, 0.058, 0.14, 0.28, 0.625, 4, 10, 40, 5.53e3, 821e3, 20e6]
nGroups = len(groups)-1

mgxs_lib = openmc.mgxs.Library(geometry)
mgxs_lib.energy_groups = openmc.mgxs.EnergyGroups(groups)
mgxs_lib.correction = None
mgxs_lib.mgxs_types = ('total','absorption','nu-fission',\
                       'fission','consistent nu-scatter matrix','chi')
mgxs_lib.domain_type = 'cell'
mgxs_lib.domains = geometry.get_all_material_cells().values()
mgxs_lib.build_library()





tallies = openmc.Tallies()
mgxs_lib.add_to_tallies_file(tallies)


flux_tally = openmc.Tally(name='flux')
energy_filter = openmc.EnergyFilter(groups)
flux_tally.filters = [openmc.CellFilter(mCells+fCells),openmc.EnergyFilter(groups)]
#flux_tally.scores = ['flux','nu-fission']
flux_tally.scores = ['flux']
tallies.append(flux_tally)


###############################################################################
###############################################################################
# Extract all Cells filled by Materials
openmc_cells = geometry.get_all_material_cells().values()

# Create dictionary to store multi-group cross sections for all cells
xs_library = {}

# Instantiate 8-group cross sections for each cell
for cell in openmc_cells:
    xs_library[cell.id] = {}
    xs_library[cell.id]['total'] = openmc.mgxs.TotalXS(groups=mgxs_lib.energy_groups)
    xs_library[cell.id]['fission'] = openmc.mgxs.FissionXS(groups=mgxs_lib.energy_groups)
    xs_library[cell.id]['nu-fission'] = openmc.mgxs.FissionXS(groups=mgxs_lib.energy_groups, nu=True)
    xs_library[cell.id]['scatter'] = openmc.mgxs.ScatterMatrixXS(groups=mgxs_lib.energy_groups)
    xs_library[cell.id]['chi'] = openmc.mgxs.Chi(groups=mgxs_lib.energy_groups)


# Iterate over all cells and cross section types
for cell in openmc_cells:
    for rxn_type in xs_library[cell.id]:

        # Set the cross sections domain to the cell
        xs_library[cell.id][rxn_type].domain = cell
        
        # Tally cross sections by nuclide
        xs_library[cell.id][rxn_type].by_nuclide = True
                
        # Add OpenMC tallies to the tallies file for XML generation
        for tally in xs_library[cell.id][rxn_type].tallies.values():
            tallies.append(tally, merge=True)
###############################################################################
###############################################################################







tallies.export_to_xml()


openmc.run()



sp = openmc.StatePoint('statepoint.100.h5')
mgxs_lib.load_from_statepoint(sp)
flux_tally=sp.get_tally(name='flux')


for i in range(9): mgxs_lib.domains[i].name   = 'fuel' + str(i)
for i in range(9): mgxs_lib.domains[9+i].name = 'mod'  + str(i)

mgxs_file = mgxs_lib.create_mg_library(xs_type='macro',xsdata_names=[       \
    'mod0', 'mod1', 'mod2', 'mod3', 'mod4', 'mod5', 'mod6', 'mod7', 'mod8', \
    'fuel0','fuel1','fuel2','fuel3','fuel4','fuel5','fuel6','fuel7','fuel8'])

mgxs_file.export_to_hdf5()

fDatas = [ mgxs_file.get_by_name('fuel'+str(i)) for i in range(9) ]
mDatas = [ mgxs_file.get_by_name('mod' +str(i)) for i in range(9) ]


# Make sure that total = scatter + absorption for both fuel and mod
for fData in fDatas:
    totalScatt = sum([fData.scatter_matrix[0][0][g][0] for g in range(nGroups)])
    assert(abs(fData.total[0][0]-(fData.absorption[0][0] + totalScatt)) < 1e-12)
for mData in mDatas:
    totalScatt = sum([mData.scatter_matrix[0][0][g][0] for g in range(nGroups)])
    assert(abs(mData.total[0][0]-(mData.absorption[0][0] + totalScatt)) < 1e-12)

###############################################################################
###############################################################################
for cell in openmc_cells:
    for rxn_type in xs_library[cell.id]:
        xs_library[cell.id][rxn_type].load_from_statepoint(sp)

nufission = xs_library[fCells[0].id]['nu-fission']
#nufission.print_xs(xs_type='macro', nuclides=['U235', 'U238','U234','U236','O16','O17'])
u234_nu_fission = nufission.get_xs(xs_type='macro', nuclides=['U234'])
u234_nu_fission_to_write = [float('%.8E'%x) for x in u234_nu_fission]
u235_nu_fission = nufission.get_xs(xs_type='macro', nuclides=['U235'])
u235_nu_fission_to_write = [float('%.8E'%x) for x in u235_nu_fission]
u236_nu_fission = nufission.get_xs(xs_type='macro', nuclides=['U236'])
u236_nu_fission_to_write = [float('%.8E'%x) for x in u236_nu_fission]
u238_nu_fission = nufission.get_xs(xs_type='macro', nuclides=['U238'])
u238_nu_fission_to_write = [float('%.8E'%x) for x in u238_nu_fission]
o16_nu_fission = nufission.get_xs(xs_type='macro', nuclides=['O16'])
o16_nu_fission_to_write = [float('%.8E'%x) for x in o16_nu_fission]
o17_nu_fission = nufission.get_xs(xs_type='macro', nuclides=['O17'])
o17_nu_fission_to_write = [float('%.8E'%x) for x in o17_nu_fission]

total= xs_library[fCells[0].id]['total']
#total.print_xs(xs_type='macro', nuclides=['U235', 'U238','U234','U236','O16','O17'])
u234_total = total.get_xs(xs_type='macro', nuclides=['U234'])
u234_total_to_write = [float('%.8E'%x) for x in u234_total]
u235_total = total.get_xs(xs_type='macro', nuclides=['U235'])
u235_total_to_write = [float('%.8E'%x) for x in u235_total]
u236_total = total.get_xs(xs_type='macro', nuclides=['U236'])
u236_total_to_write = [float('%.8E'%x) for x in u236_total]
u238_total = total.get_xs(xs_type='macro', nuclides=['U238'])
u238_total_to_write = [float('%.8E'%x) for x in u238_total]
o16_total = total.get_xs(xs_type='macro', nuclides=['O16'])
o16_total_to_write = [float('%.8E'%x) for x in o16_total]
o17_total = total.get_xs(xs_type='macro', nuclides=['O17'])
o17_total_to_write = [float('%.8E'%x) for x in o17_total]


scatter = xs_library[fCells[0].id]['scatter']
#scatter.print_xs(xs_type='macro', nuclides=['U235', 'U238','U234','U236','O16','O17'])
u234_scatter = scatter.get_xs(xs_type='macro', nuclides=['U234'])
u234_scatter_to_write = [[float('%.8E'%x) for x in y] for y in u234_scatter]
u235_scatter = scatter.get_xs(xs_type='macro', nuclides=['U235'])
u235_scatter_to_write = [[float('%.8E'%x) for x in y] for y in u235_scatter]
u236_scatter = scatter.get_xs(xs_type='macro', nuclides=['U236'])
u236_scatter_to_write = [[float('%.8E'%x) for x in y] for y in u236_scatter]
u238_scatter = scatter.get_xs(xs_type='macro', nuclides=['U238'])
u238_scatter_to_write = [[float('%.8E'%x) for x in y] for y in u238_scatter]
o16_scatter = scatter.get_xs(xs_type='macro', nuclides=['O16'])
o16_scatter_to_write = [[float('%.8E'%x) for x in y] for y in o16_scatter]
o17_scatter = scatter.get_xs(xs_type='macro', nuclides=['O17'])
o17_scatter_to_write = [[float('%.8E'%x) for x in y] for y in o17_scatter]

chi= xs_library[fCells[0].id]['chi']
#chi.print_xs(xs_type='macro', nuclides=['U235', 'U238','U234','U236','O16','O17'])
u234_chi = chi.get_xs(xs_type='macro', nuclides=['U234'])
u234_chi_to_write = [float('%.8E'%x) for x in u234_chi]
u235_chi = chi.get_xs(xs_type='macro', nuclides=['U235'])
u235_chi_to_write = [float('%.8E'%x) for x in u235_chi]
u236_chi = chi.get_xs(xs_type='macro', nuclides=['U236'])
u236_chi_to_write = [float('%.8E'%x) for x in u236_chi]
u238_chi = chi.get_xs(xs_type='macro', nuclides=['U238'])
u238_chi_to_write = [float('%.8E'%x) for x in u238_chi]
o16_chi = chi.get_xs(xs_type='macro', nuclides=['O16'])
o16_chi_to_write = [float('%.8E'%x) for x in o16_chi]
o17_chi = chi.get_xs(xs_type='macro', nuclides=['O17'])
o17_chi_to_write = [float('%.8E'%x) for x in o17_chi]





###############################################################################
###############################################################################



f = open("XS_nuclideSpecific.py","w+")
#
for i in range(9):
    f.write("u234_nuFission = "+str(u234_nu_fission_to_write)+"\n")
    f.write("u235_nuFission = "+str(u235_nu_fission_to_write)+"\n")
    f.write("u236_nuFission = "+str(u236_nu_fission_to_write)+"\n")
    f.write("u238_nuFission = "+str(u238_nu_fission_to_write)+"\n")
    f.write("o16_nuFission = "+str(o16_nu_fission_to_write)+"\n")
    f.write("o17_nuFission = "+str(o17_nu_fission_to_write)+"\n")
    f.write("\n\n")
    f.write("u234_total = "+str(u234_total_to_write)+"\n")
    f.write("u235_total = "+str(u235_total_to_write)+"\n")
    f.write("u236_total = "+str(u236_total_to_write)+"\n")
    f.write("u238_total = "+str(u238_total_to_write)+"\n")
    f.write("o16_total = "+str(o16_total_to_write)+"\n")
    f.write("o17_total = "+str(o17_total_to_write)+"\n")
    f.write("\n\n")
    f.write("u234_scatter = "+str(u234_scatter_to_write)+"\n")
    f.write("u235_scatter = "+str(u235_scatter_to_write)+"\n")
    f.write("u236_scatter = "+str(u236_scatter_to_write)+"\n")
    f.write("u238_scatter = "+str(u238_scatter_to_write)+"\n")
    f.write("o16_scatter = "+str(o16_scatter_to_write)+"\n")
    f.write("o17_scatter = "+str(o17_scatter_to_write)+"\n")
    f.write("\n\n")
    f.write("u234_chi = "+str(u234_chi_to_write)+"\n")
    f.write("u235_chi = "+str(u235_chi_to_write)+"\n")
    f.write("u236_chi = "+str(u236_chi_to_write)+"\n")
    f.write("u238_chi = "+str(u238_chi_to_write)+"\n")
    f.write("o16_chi = "+str(o16_chi_to_write)+"\n")
    f.write("o17_chi = "+str(o17_chi_to_write)+"\n")
    f.write("\n\n")

    f.write("u234_all = [u234_nuFission,u234_total,u234_scatter,u234_chi]\n")
    f.write("u235_all = [u235_nuFission,u235_total,u235_scatter,u235_chi]\n")
    f.write("u236_all = [u236_nuFission,u236_total,u236_scatter,u236_chi]\n")
    f.write("u238_all = [u238_nuFission,u238_total,u238_scatter,u238_chi]\n")
    f.write("o16_all = [o16_nuFission,o16_total,o16_scatter,o16_chi]\n")
    f.write("o17_all = [o17_nuFission,o17_total,o17_scatter,o17_chi]\n")

    break
    f.write("\n\n")



f.close()



