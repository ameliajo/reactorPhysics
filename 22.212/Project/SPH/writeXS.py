

def writeXS( sph, nGroups, nPins, filename,                                  \
             fuelTotal, fuelAbsorption, fuelNuFission, fuelChi, fuelScatter, \
             modTotal,  modAbsorption,  modNuFission,  modChi,  modScatter ):
    f = open(filename,"w+")

    for g in range(nGroups): sph[0][g] = 1.0 if sph[0][g] < 1e-10 else sph[0][g]
    for g in range(nGroups): sph[1][g] = 1.0 if sph[1][g] < 1e-10 else sph[1][g]

    f_total_new = [fuelTotal[g]/sph[0][g] for g in range(nGroups)]
    f_absorption_new = [fuelAbsorption[g]/sph[0][g] for g in range(nGroups)]
    f_nuFission_new = [fuelNuFission[g]/sph[0][g] for g in range(nGroups)]
    #f_chi_new = [fuelChi[g]/sph[0][g] for g in range(nGroups)]
    f_chi_new = fuelChi
    f_scatter_new = [[fuelScatter[g][gp]/sph[0][g] for gp in range(nGroups)] for g in range(nGroups)]
    
    i = 0
    f.write("fuelTotal"+str(i)+"  = "+str([float("%.8f"%x) for x in f_total_new])+"\n")
    f.write("fuelAbsorption"+str(i)+"  = "+str([float("%.8f"%x) for x in f_absorption_new])+"\n")
    f.write("fuelNuFission"+str(i)+"  = "+str([float("%.8f"%x) for x in f_nuFission_new])+"\n")
    f.write("fuelChi"+str(i)+"  = "+str([float("%.8f"%x) for x in f_chi_new])+"\n")
    f.write("fuelScatter"+str(i)+" = "+str([[float("%.8f"%f_scatter_new[g][gp]) for gp in range(nGroups)] for g in range(nGroups)])+"\n")
    
    """
    m_total_new = [modTotal[g] /sph[1][g] for g in range(nGroups)]
    m_absorption_new = [modAbsorption[g]/sph[1][g] for g in range(nGroups)]
    m_nuFission_new = [modNuFission[g]/sph[1][g] for g in range(nGroups)]
    #m_chi_new = [modChi[g]/sph[1][g] for g in range(nGroups)]
    m_chi_new = modChi
    m_scatter_new = [[modScatter[g][gp]/sph[1][g] for gp in range(nGroups)] for g in range(nGroups)]
    """

    m_total_new = modTotal
    m_absorption_new = modAbsorption
    m_nuFission_new = modNuFission
    m_chi_new = modChi
    m_scatter_new = modScatter
 
        
    f.write("modTotal"+str(i)+"  = "+str([float("%.8f"%x) for x in m_total_new])+"\n")
    f.write("modAbsorption"+str(i)+"  = "+str([float("%.8f"%x) for x in m_absorption_new])+"\n")
    f.write("modNuFission"+str(i)+"  = "+str([float("%.8f"%x) for x in m_nuFission_new])+"\n")
    f.write("modChi"+str(i)+"  = "+str([float("%.8f"%x) for x in m_chi_new])+"\n")
    f.write("modScatter"+str(i)+" = "+str([[float("%.8f"%m_scatter_new[g][gp]) for gp in range(nGroups)] for g in range(nGroups)])+"\n")
 

    f.close()


    """
    for i in range(nPins):

        f_total_new = [fuelTotal[i][g] /sph[i][0][g] for g in range(nGroups)]
        f_absorption_new = [fuelAbsorption[i][g]/sph[i][0][g] for g in range(nGroups)]
        f_nuFission_new = [fuelNuFission[i][g]/sph[i][0][g] for g in range(nGroups)]
        #f_chi_new = [fuelChi[i][g]/sph[i][0][g] for g in range(nGroups)]
        f_chi_new = fuelChi[i]
        f_scatter_new = [[fuelScatter[i][g][gp]/sph[i][0][g] for gp in range(nGroups)] for g in range(nGroups)]
    
        f.write("fuelTotal"+str(i)+"  = "+str([float("%.8f"%x) for x in f_total_new])+"\n")
        f.write("fuelAbsorption"+str(i)+"  = "+str([float("%.8f"%x) for x in f_absorption_new])+"\n")
        f.write("fuelNuFission"+str(i)+"  = "+str([float("%.8f"%x) for x in f_nuFission_new])+"\n")
        f.write("fuelChi"+str(i)+"  = "+str([float("%.8f"%x) for x in f_chi_new])+"\n")
        f.write("fuelScatter"+str(i)+" = "+str([[float("%.8f"%f_scatter_new[g][gp]) for gp in range(nGroups)] for g in range(nGroups)])+"\n")
    

        m_total_new = [modTotal[i][g] /sph[i][1][g] for g in range(nGroups)]
        m_absorption_new = [modAbsorption[i][g]/sph[i][1][g] for g in range(nGroups)]
        m_nuFission_new = [modNuFission[i][g]/sph[i][1][g] for g in range(nGroups)]
        #m_chi_new = [modChi[i][g]/sph[i][1][g] for g in range(nGroups)]
        m_chi_new = modChi[i]
        m_scatter_new = [[modScatter[i][g][gp]/sph[i][1][g] for gp in range(nGroups)] for g in range(nGroups)]
        
        f.write("modTotal"+str(i)+"  = "+str([float("%.8f"%x) for x in m_total_new])+"\n")
        f.write("modAbsorption"+str(i)+"  = "+str([float("%.8f"%x) for x in m_absorption_new])+"\n")
        f.write("modNuFission"+str(i)+"  = "+str([float("%.8f"%x) for x in m_nuFission_new])+"\n")
        f.write("modChi"+str(i)+"  = "+str([float("%.8f"%x) for x in m_chi_new])+"\n")
        f.write("modScatter"+str(i)+" = "+str([[float("%.8f"%m_scatter_new[g][gp]) for gp in range(nGroups)] for g in range(nGroups)])+"\n")
 

    f.close()

    """


