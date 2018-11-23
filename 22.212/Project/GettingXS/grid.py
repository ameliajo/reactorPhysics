

#get_ipython().magic('matplotlib inline')
import openmc
import openmc.mgxs as mgxs
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


uo2 = openmc.Material(1,"uo2")
uo2.add_element('U', 1.0, enrichment=3.0)
uo2.add_element('O', 2.0)
uo2.set_density('g/cc', 10.0)

uo2_b = openmc.Material(2,"uo2")
uo2_b.add_element('U', 1.0, enrichment=3.0)
uo2_b.add_element('O', 2.0)
uo2_b.set_density('g/cc', 10.0)


water = openmc.Material(3, "h2o")
water.add_element('H', 2.0)
water.add_element('O', 1.0)
water.set_density('g/cm3', 1.0)

water_b = openmc.Material(4, "h2o")
water_b.add_element('H', 2.0)
water_b.add_element('O', 1.0)
water_b.set_density('g/cm3', 1.0)


mats = openmc.Materials([uo2, uo2_b, water, water_b])

mats.export_to_xml()


##################################################################
# GEOMETRY
##################################################################
L = 1.26
pitch = L*3

fuel_or_vec = [ openmc.ZCylinder(R=0.39, x0=0.5*L+i*L, y0=0.5*L+j*L) for j in range(3) for i in range(3) ]


F_cells = [ openmc.Cell(i,'fuel'+str(i),fill=uo2,region=-fuel_or_vec[i]) for i in range(9)]

x1 = openmc.XPlane(x0=0.0, boundary_type='reflective')
x2 = openmc.XPlane(x0=1*L)
x3 = openmc.XPlane(x0=2*L)
x4 = openmc.XPlane(x0=3*L, boundary_type='reflective')

y1 = openmc.YPlane(y0=0.0, boundary_type='reflective')
y2 = openmc.YPlane(y0=1*L)
y3 = openmc.YPlane(y0=2*L)
y4 = openmc.YPlane(y0=3*L, boundary_type='reflective')

waterRegion = +x1 & -x4 & +y1 & -y4 &                                 \
              +fuel_or_vec[0] &  +fuel_or_vec[1] &  +fuel_or_vec[2] & \
              +fuel_or_vec[3] &  +fuel_or_vec[4] &  +fuel_or_vec[5] & \
              +fuel_or_vec[6] &  +fuel_or_vec[7] &  +fuel_or_vec[8] 

xP = [x1,x2,x3,x4]
yP = [y1,y2,y3,y4]

M_cells = []
count = 0
for j in range(3):
    for i in range(3):
        mRegion = waterRegion & +xP[i] & -xP[i+1] & +yP[j] & -yP[j+1]
        M_cells.append(openmc.Cell(70+count,'mod'+str(i),fill=water,region=mRegion))
        count += 1

#special = 7
#F_cells[special].fill=uo2_b
#M_cells[special].fill=water_b

root = openmc.Universe(cells=(                                              \
    M_cells[0], M_cells[1], M_cells[2], M_cells[3], M_cells[4], M_cells[5], \
    M_cells[6], M_cells[7], M_cells[8], F_cells[0], F_cells[1], F_cells[2], \
    F_cells[3], F_cells[4], F_cells[5], F_cells[6], F_cells[7], F_cells[8] ))


geom = openmc.Geometry(root)
geom.export_to_xml()


##################################################################
# SETTINGS
##################################################################

settings = openmc.Settings()
#settings.batches = 500
#ettings.inactive = 50
#ettings.particles = 10000
settings.batches = 200
settings.inactive = 10
settings.particles = 1000

settings.output = {'tallies': True}

bounds = [0.0, 0.0, 0.0, pitch, pitch, pitch]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:])
settings.source = openmc.source.Source(space=uniform_dist)

settings.export_to_xml()



# Instantiate a 2-group EnergyGroups object
groups = mgxs.EnergyGroups()
#groups.group_edges = np.array([0., 0.625, 20.0e6])
groups.group_edges = np.array([0.0, 0.058, 0.14, 0.28, 0.625, 4, 10, 40, 5.53e3, 821e3, 20e6])

totalFuel = [mgxs.TotalXS(        domain=F_cell, groups=groups) for F_cell in F_cells]
absorFuel = [mgxs.AbsorptionXS(   domain=F_cell, groups=groups) for F_cell in F_cells]
scattFuel = [mgxs.ScatterMatrixXS(domain=F_cell, groups=groups) for F_cell in F_cells]
fizzzFuel = [mgxs.FissionXS(      domain=F_cell, groups=groups, nu=True) for F_cell in F_cells]
chiiiFuel = [mgxs.Chi(            domain=F_cell, groups=groups) for F_cell in F_cells]

totalMod  = [mgxs.TotalXS(        domain=M_cell, groups=groups) for M_cell in M_cells]
absorMod  = [mgxs.AbsorptionXS(   domain=M_cell, groups=groups) for M_cell in M_cells]
scattMod  = [mgxs.ScatterMatrixXS(domain=M_cell, groups=groups) for M_cell in M_cells]


tallies_file = openmc.Tallies()

for totalF in totalFuel: tallies_file += totalF.tallies.values()
for absorF in absorFuel: tallies_file += absorF.tallies.values()
for scattF in scattFuel: tallies_file += scattF.tallies.values()
for fizzzF in fizzzFuel: tallies_file += fizzzF.tallies.values()
for chiiiF in chiiiFuel: tallies_file += chiiiF.tallies.values()

for totalM in totalMod: tallies_file += totalM.tallies.values()
for absorM in absorMod: tallies_file += absorM.tallies.values()
for scattM in scattMod: tallies_file += scattM.tallies.values()


tallies_file.export_to_xml()

openmc.run()

sp = openmc.StatePoint('statepoint.200.h5')

for i in range(len(totalFuel)): totalFuel[i].load_from_statepoint(sp)
for i in range(len(absorFuel)): absorFuel[i].load_from_statepoint(sp)
for i in range(len(scattFuel)): scattFuel[i].load_from_statepoint(sp)
for i in range(len(fizzzFuel)): fizzzFuel[i].load_from_statepoint(sp)
for i in range(len(chiiiFuel)): chiiiFuel[i].load_from_statepoint(sp)

for i in range(len(totalMod)): totalMod[i].load_from_statepoint(sp)
for i in range(len(absorMod)): absorMod[i].load_from_statepoint(sp)
for i in range(len(scattMod)): scattMod[i].load_from_statepoint(sp)


flux = [[float(str("%0.6f"%y[0][0])) for y in totalFuel[i].tallies['flux'].get_slice(scores=['flux']).mean] for i in range(9)]
mod  = [[float(str("%0.6f"%y[0][0])) for y in totalMod[i].tallies['flux'].get_slice(scores=['flux']).mean] for i in range(9)]

print(flux)
print(mod)

f = open("fuelXS.txt","w+")
f_py = open("fuelXS.py","w+")

for i in range(9):
    df = totalFuel[i].get_pandas_dataframe()
    coeff = F_cells[i].region.surface.coefficients
    f.write("FUEL CELL "+str(i+1)+", with center at ("+str(coeff['x0'])+","+\
            str(coeff['y0'])+")\n")
    f.write("-------------------------------------------------------------------------------\n\n\n")
    f.write("TOTAL FUEL CELL "+str(i+1)+", with center at ("+ \
            str(coeff['x0'])+","+str(coeff['y0'])+")\n")
    f.write(str(df))
    f.write("\n\n")
    f_py.write("totalFuel"+str(i)+" = "+str([float("%.5f" % df) for df in df['mean'].values.tolist()])+"\n")

    df = absorFuel[i].get_pandas_dataframe()
    coeff = F_cells[i].region.surface.coefficients
    f.write("ABSOR FUEL CELL "+str(i+1)+", with center at ("+ \
            str(coeff['x0'])+","+str(coeff['y0'])+")\n")
    f.write(str(df))
    f.write("\n\n")
    f_py.write("absorFuel"+str(i)+" = "+str([float("%.5f" % df) for df in df['mean'].values.tolist()])+"\n")

    df = scattFuel[i].get_pandas_dataframe()
    coeff = F_cells[i].region.surface.coefficients

    f.write("SCATTER FUEL CELL "+str(i+1)+", with center at ("+ \
            str(coeff['x0'])+","+str(coeff['y0'])+")\n")
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        f.write(str(df))
    f.write("\n\n")
    f_py.write("scattFuel"+str(i)+" = "+str([float("%.5f" % df) for df in df['mean'].values.tolist()])+"\n")


    df = fizzzFuel[i].get_pandas_dataframe()
    coeff = F_cells[i].region.surface.coefficients
    f.write("NU-FISSION FUEL CELL "+str(i+1)+", with center at ("+ \
            str(coeff['x0'])+","+str(coeff['y0'])+")\n")
    f.write(str(df))
    f_py.write("fizzzFuel"+str(i)+" = "+str([float("%.5f" % df) for df in df['mean'].values.tolist()])+"\n")
    f.write("\n\n\n\n\n")

    df = chiiiFuel[i].get_pandas_dataframe()
    coeff = F_cells[i].region.surface.coefficients
    f.write("CHI FUEL CELL "+str(i+1)+", with center at ("+ \
            str(coeff['x0'])+","+str(coeff['y0'])+")\n")
    f.write(str(df))
    f_py.write("chiiiFuel"+str(i)+" = "+str([float("%.5f" % df) for df in df['mean'].values.tolist()])+"\n")
    f.write("\n\n\n\n\n")






f.close()
f_py.close()



f= open("modXS.txt","w+")
f_py= open("modXS.py","w+")

for i in range(9):
    df = totalMod[i].get_pandas_dataframe()
    coeff = F_cells[i].region.surface.coefficients
    f.write("FUEL CELL "+str(i+1)+", with center at ("+str(coeff['x0'])+","+\
            str(coeff['y0'])+")\n")
    f.write("-------------------------------------------------------------------------------\n\n\n")
    f.write("TOTAL FUEL CELL "+str(i+1)+", with center at ("+ \
            str(coeff['x0'])+","+str(coeff['y0'])+")\n")
    f.write(str(df))
    f.write("\n\n")
    f_py.write("totalMod"+str(i)+" = "+str([float("%.5f" % df) for df in df['mean'].values.tolist()])+"\n")

    df = absorMod[i].get_pandas_dataframe()
    coeff = F_cells[i].region.surface.coefficients
    f.write("ABSOR FUEL CELL "+str(i+1)+", with center at ("+ \
            str(coeff['x0'])+","+str(coeff['y0'])+")\n")
    f.write(str(df))
    f.write("\n\n")
    f_py.write("absorMod"+str(i)+" = "+str([float("%.5f" % df) for df in df['mean'].values.tolist()])+"\n")

    df = scattMod[i].get_pandas_dataframe()
    coeff = F_cells[i].region.surface.coefficients

    f.write("SCATTER FUEL CELL "+str(i+1)+", with center at ("+ \
            str(coeff['x0'])+","+str(coeff['y0'])+")\n")
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        f.write(str(df))
    f.write("\n\n")
    f_py.write("scattMod"+str(i)+" = "+str([float("%.5f" % df) for df in df['mean'].values.tolist()])+"\n")


    f.write("\n\n")
    f.write("\n\n")
    f.write("\n\n")




f.close()
f_py.close()



##################################################################
# PLOT
##################################################################
import sys
if (len(sys.argv) > 1):
    if (sys.argv[1] == 'plot'):

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

