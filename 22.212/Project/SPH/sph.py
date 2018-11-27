import subprocess
import matplotlib.pyplot as plt
import sys

if len(sys.argv) > 1:
    if sys.argv[1] == "full":
        subprocess.run(['ls'])
        subprocess.run(['python','../GettingXS/grid_3x3.py'])
        subprocess.run(['cp','XS.py','../MOC'])
        subprocess.run(['rm','geometry.xml','materials.xml','mgxs.h5','settings.xml','statepoint.100.h5','summary.h5','tallies.out','tallies.xml'])
        subprocess.run(['ls'])
        subprocess.run(['python','../MOC/moc.py'])


from fluxMOC import *
from fluxMC  import *
from XS import *


MOC_modFlux = [ MOC_modFlux0, MOC_modFlux1, MOC_modFlux2, MOC_modFlux3, 
                MOC_modFlux4, MOC_modFlux5, MOC_modFlux6, MOC_modFlux7, MOC_modFlux8 ]
MC_modFlux  = [ MC_modFlux0, MC_modFlux1, MC_modFlux2, MC_modFlux3, MC_modFlux4, 
                MC_modFlux5, MC_modFlux6, MC_modFlux7, MC_modFlux8 ]

MOC_fuelFlux = [ MOC_fuelFlux0, MOC_fuelFlux1, MOC_fuelFlux2, MOC_fuelFlux3, 
                 MOC_fuelFlux4, MOC_fuelFlux5, MOC_fuelFlux6, MOC_fuelFlux7, MOC_fuelFlux8 ]
MC_fuelFlux  = [ MC_fuelFlux0, MC_fuelFlux1, MC_fuelFlux2, MC_fuelFlux3, 
                 MC_fuelFlux4, MC_fuelFlux5, MC_fuelFlux6, MC_fuelFlux7, MC_fuelFlux8 ]

MC_sum  = 0.0
MOC_sum = 0.0
for i in range(9):
    MC_sum  += sum(MC_modFlux[i])  + sum(MC_fuelFlux[i])
    MOC_sum += sum(MOC_modFlux[i]) + sum(MOC_fuelFlux[i])


print([MOC_fuelFlux0[i]/MC_fuelFlux0[i] for i in range(len(MOC_fuelFlux0))])
print([MOC_fuelFlux1[i]/MC_fuelFlux0[i] for i in range(len(MOC_fuelFlux1))])
print([MOC_fuelFlux2[i]/MC_fuelFlux0[i] for i in range(len(MOC_fuelFlux2))])
print([MOC_fuelFlux3[i]/MC_fuelFlux0[i] for i in range(len(MOC_fuelFlux3))])

nGroups = len(MC_modFlux0)
invSumMC  = 1.0/MC_sum
invSumMOC = 1.0/MOC_sum


for i in range(9):

    MOC_mod = [x * invSumMOC for x in MOC_modFlux[i]]
    MC_mod  = [x * invSumMC  for x in MC_modFlux[i]]

    MOC_fuel = [x * invSumMOC for x in MOC_fuelFlux[i]]
    MC_fuel  = [x * invSumMC  for x in MC_fuelFlux[i]]
    
    plt.plot(MOC_mod,'ro-',label='MOC moderator'+str(i))
    plt.plot(MC_mod, 'rx-',label='MC  moderator'+str(i))
    plt.plot(MOC_fuel,'bo-',label='MOC fuel'+str(i))
    plt.plot(MC_fuel, 'bx-',label='MC  fuel'+str(i))


    break


sph = []


for cellNum in range(9):
    sph.append([[],[]])
    for g in range(nGroups):
        sph_mod  = MC_modFlux[cellNum][g] / MOC_modFlux[cellNum][g]
        sph_fuel = MC_modFlux[cellNum][g] / MOC_modFlux[cellNum][g]
        sph[cellNum][0].append(sph_mod)
        sph[cellNum][1].append(sph_fuel)


plt.legend(loc='best')



fuelTotal = [ fuelTotal0, fuelTotal1, fuelTotal2, fuelTotal3, fuelTotal4, fuelTotal5, fuelTotal6, fuelTotal7, fuelTotal8]
fuelAbsorption = [ fuelAbsorption0, fuelAbsorption1, fuelAbsorption2, fuelAbsorption3, fuelAbsorption4, fuelAbsorption5, fuelAbsorption6, fuelAbsorption7, fuelAbsorption8]
fuelChi = [ fuelChi0, fuelChi1, fuelChi2, fuelChi3, fuelChi4, fuelChi5, fuelChi6, fuelChi7, fuelChi8]
fuelNuFission = [ fuelNuFission0, fuelNuFission1, fuelNuFission2, fuelNuFission3, fuelNuFission4, fuelNuFission5, fuelNuFission6, fuelNuFission7, fuelNuFission8]
fuelScatter = [ fuelScatter0, fuelScatter1, fuelScatter2, fuelScatter3, fuelScatter4, fuelScatter5, fuelScatter6, fuelScatter7, fuelScatter8]


modTotal = [ modTotal0, modTotal1, modTotal2, modTotal3, modTotal4, modTotal5, modTotal6, modTotal7, modTotal8]
modAbsorption = [ modAbsorption0, modAbsorption1, modAbsorption2, modAbsorption3, modAbsorption4, modAbsorption5, modAbsorption6, modAbsorption7, modAbsorption8]
modChi = [ modChi0, modChi1, modChi2, modChi3, modChi4, modChi5, modChi6, modChi7, modChi8]
modNuFission = [ modNuFission0, modNuFission1, modNuFission2, modNuFission3, modNuFission4, modNuFission5, modNuFission6, modNuFission7, modNuFission8]
modScatter = [ modScatter0, modScatter1, modScatter2, modScatter3, modScatter4, modScatter5, modScatter6, modScatter7, modScatter8]




f = open("XS2.py","w+")
for i in range(9):

    f_total_new = [fuelTotal[i][g] /sph[i][0][g] for g in range(nGroups)]
    f_absorption_new = [fuelAbsorption[i][g]/sph[i][0][g] for g in range(nGroups)]
    f_nuFission_new = [fuelNuFission[i][g]/sph[i][0][g] for g in range(nGroups)]
    f_chi_new = [fuelChi[i][g]/sph[i][0][g] for g in range(nGroups)]
    f_scatter_new = [[fuelScatter[i][g][gp]/sph[i][0][g] for gp in range(nGroups)] for g in range(nGroups)]
    
    f.write("fuelTotal"+str(i)+"  = "+str([float("%.8f"%x) for x in f_total_new])+"\n")
    f.write("fuelAbsorption"+str(i)+"  = "+str([float("%.8f"%x) for x in f_absorption_new])+"\n")
    f.write("fuelNuFission"+str(i)+"  = "+str([float("%.8f"%x) for x in f_nuFission_new])+"\n")
    f.write("fuelChi"+str(i)+"  = "+str([float("%.8f"%x) for x in f_chi_new])+"\n")
    f.write("fuelScatter"+str(i)+" = "+str([[float("%.8f"%f_scatter_new[g][gp]) for gp in range(nGroups)] for g in range(nGroups)])+"\n")
    

    m_total_new = [modTotal[i][g] /sph[i][1][g] for g in range(nGroups)]
    m_absorption_new = [modAbsorption[i][g]/sph[i][1][g] for g in range(nGroups)]
    m_nuFission_new = [modNuFission[i][g]/sph[i][1][g] for g in range(nGroups)]
    m_chi_new = [modChi[i][g]/sph[i][1][g] for g in range(nGroups)]
    m_scatter_new = [[modScatter[i][g][gp]/sph[i][1][g] for gp in range(nGroups)] for g in range(nGroups)]
    
    f.write("modTotal"+str(i)+"  = "+str([float("%.8f"%x) for x in m_total_new])+"\n")
    f.write("modAbsorption"+str(i)+"  = "+str([float("%.8f"%x) for x in m_absorption_new])+"\n")
    f.write("modNuFission"+str(i)+"  = "+str([float("%.8f"%x) for x in m_nuFission_new])+"\n")
    f.write("modChi"+str(i)+"  = "+str([float("%.8f"%x) for x in m_chi_new])+"\n")
    f.write("modScatter"+str(i)+" = "+str([[float("%.8f"%m_scatter_new[g][gp]) for gp in range(nGroups)] for g in range(nGroups)])+"\n")
 
    print("SPH FACTOR",sph[0][0][0])
    print("SPH FACTOR",sph[1][1][1])
    print("SPH FACTOR",sph)
    break


f.close()



"""
for i in range(9):
    f.write("fuelTotal"+str(i)+"  = "+str([float("%.8f"%flux) for flux in phi[i][0]])+"\n")
    f.write("modTotal"+str(i)+" = "+str([float("%.8f"%flux) for flux in phi[i][1]])+"\n")
    f.write("\n\n")


"""





















