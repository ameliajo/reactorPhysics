from math import pi
import openmc
import openmc.mgxs
import openmc.model
import numpy as np


radius_fuel = 0.39128
pitch = 1.26

# Basic materials
uo2 = openmc.Material(name='fuel')
uo2.add_element('U', 1, enrichment=3.2)
uo2.add_element('O', 2)
uo2.add_element('Gd', 0.0007)
uo2.set_density('g/cc', 10.341)


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
settings.particles = 500



space = openmc.stats.Box((0.0, 0.0, 0.0),(3.0*pitch, 3.0*pitch, 0))
settings.source = openmc.Source(space=space)
settings.export_to_xml()

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
flux_tally.filters = [openmc.CellFilter(mCells+fCells)]
flux_tally.filters.append(energy_filter)
#flux_tally.scores = ['flux','nu-fission']
flux_tally.scores = ['flux']
tallies.append(flux_tally)
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
        openmc.plot_geometry()





f = open("XS.py","w+")

for i in range(9):
    fData = fDatas[i]
    mData = mDatas[i]
    coeff = fCells[i].region.surface.coefficients
    x0 = str("%.5f" % coeff['x0'])
    y0 = str("%.5f" % coeff['y0'])

    f.write("# CELL "+str(i+1)+", with center at ("+x0+","+y0+")\n")
    f.write("# -----------------------------------------------------------------------------\n\n")

    f.write("fuelTotal"+str(i)+" = "+str([float("%.8f"%f) for f in fData.total[0]])+"\n")
    f.write("fuelAbsorption"+str(i)+" = "+str([float("%.8f"%f) for f in fData.absorption[0]])+"\n")
    f.write("fuelNuFission"+str(i)+" = "+str([float("%.8f"%f) for f in fData.nu_fission[0]])+"\n")
    f.write("fuelChi"+str(i)+" = "+str([float("%.8f"%f) for f in fData.chi[0]])+"\n")
    f.write("fuelScatter"+str(i)+" = "+str([[float("%.8f"%fData.scatter_matrix[0][g][gp][0]) for gp in range(nGroups)] for g in range(nGroups)])+"\n")
    f.write("# SigS[g][g'] = Scattering g->g'. So "+str(float("%.8f"%fData.scatter_matrix[0][1-1][2-1][0]))+" is scattering from 1->2"+"\n")

    f.write("\n")

    f.write("modTotal"+str(i)+" = "+str([float("%.8f"%m) for m in mData.total[0]])+"\n")
    f.write("modAbsorption"+str(i)+" = "+str([float("%.8f"%m) for m in mData.absorption[0]])+"\n")
    f.write("modNuFission"+str(i)+" = "+str([float("%.8f"%m) for m in mData.nu_fission[0]])+"\n")
    f.write("modChi"+str(i)+" = "+str([float("%.8f"%m) for m in mData.chi[0]])+"\n")
    f.write("modScatter"+str(i)+" = "+str([[float("%.8f"%mData.scatter_matrix[0][g][gp][0]) for gp in range(nGroups)] for g in range(nGroups)])+"\n")
    f.write("# SigS[g][g'] = Scattering g->g'. So "+str(float("%.8f"%mData.scatter_matrix[0][1-1][2-1][0]))+" is scattering from 1->2"+"\n")

    f.write("\n\n")



f.close()





f = open("fluxMC.py","w+")

for i in range(9):
    modFlux  = flux_tally.get_slice(filters=[openmc.CellFilter], filter_bins=[((mCells[i]).id,)])
    fuelFlux = flux_tally.get_slice(filters=[openmc.CellFilter], filter_bins=[((fCells[i]).id,)])
    f.write("MC_modFlux"+str(i)+"  = "+str([float("%.8f"%flux[0][0]) for flux in modFlux.mean])+"\n")
    f.write("MC_fuelFlux"+str(i)+" = "+str([float("%.8f"%flux[0][0]) for flux in fuelFlux.mean])+"\n")
    f.write("\n\n")
f.close()



























