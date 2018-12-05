import time
from numpy.random import random_sample as rand
from plotting import *
from ray import *
from physics import *
from materials import *

np.random.seed(42)

def runRays(n_rays, surfaces, regions, sideLen, ngroup, plot=False, 
            maxRayDist=100.0, deadzone=10, verbose=True, sphIter=0):

    start = time.perf_counter()

    ##########################################################################
    # DRAW ALL RAYS
    ##########################################################################
    rays = []
    allActiveDists = 0.0
    rays = [drawRay(Ray(sideLen,maxRayDist), surfaces, regions, deadzone) for i in range(n_rays)]
    allActiveDists = sum([ray.distActive for ray in rays])

    for region in regions:
        region.vol = region.activeDist/allActiveDists
    


    ##########################################################################
    # ITERATE  
    ##########################################################################
    oldFissSrc, k = updateQ(regions, ngroup, 1, justStarting=True)
    kVals = [k]
    converged, counter = False, 0

    while not converged:
        counter += 1
        if verbose: print('Iterations # ', counter, ' k = ', k)

        for reg in regions: reg.phi = np.zeros([ngroup])

        updatePhi(rays, ngroup, regions)

        for reg in regions:
            reg.phi = reg.phi/(reg.vol*reg.mat.SigT*allActiveDists) + reg.q/reg.mat.SigT

        newFissSrc, k = updateQ(regions, ngroup, k, oldFissSrc)
        diffFissSrc = abs(newFissSrc-oldFissSrc)/oldFissSrc
        diffK       = abs((kVals[-1]-k)/kVals[-1])
        oldFissSrc = newFissSrc

        converged = (diffK < 1e-5 and diffFissSrc < 1e-7) or (counter > 100)
        kVals.append(k)
        
    if verbose: print('k = ', k, ' after ', counter, 'iterations')

    end = time.perf_counter()
    elapsed_time = end - start

    segments = sum([len(ray.segments) for ray in rays])

    if verbose: print('Elapsed time:               ', elapsed_time)
    if verbose: print('Time per segment per group: ', elapsed_time/(segments*ngroup))





    ##########################################################################
    # WRITING FLUX TO FILE
    ##########################################################################

    #print("Writing to fluxMOC"+str(sphIter)+".py")
    f = open(str("fluxMOC"+str(sphIter)+".py"),"w+")
    f.write("MOC_k = "+str(k)+"\n")
    for i in range(9):
        f.write("MOC_modFlux"+str(i)+" = "+str([float("%.8f"%flux) for flux in regions[9+i].phi])+"\n")
        f.write("MOC_fuelFlux"+str(i)+" = "+str([float("%.8f"%flux) for flux in regions[i].phi])+"\n")
    f.write("\n\n")
    f.write("MOC_modFlux = [ MOC_modFlux0, MOC_modFlux1, MOC_modFlux2, MOC_modFlux3, MOC_modFlux4, MOC_modFlux5, MOC_modFlux6, MOC_modFlux7, MOC_modFlux8 ]\n")
    f.write("MOC_fuelFlux = [ MOC_fuelFlux0, MOC_fuelFlux1, MOC_fuelFlux2, MOC_fuelFlux3, MOC_fuelFlux4, MOC_fuelFlux5, MOC_fuelFlux6, MOC_fuelFlux7, MOC_fuelFlux8 ]")
    f.close()

    """
    for i,region in enumerate(regions):
        f.write("MOC_modFlux"+str(i)+"  = "+str([float("%.8f"%flux) for flux in region.phi])+"\n")
        f.write("MOC_fuelFlux"+str(i)+" = "+str([float("%.8f"%flux) for flux in region.phi])+"\n")
        f.write("\n\n")
    f.close()
    """





    ##########################################################################
    # PLOTTING
    ##########################################################################
    if plot:
        ktitle ='k = '+str(k)+' Rays ='+str(n_rays)
        print('Plotting tracks')
        plot_from_rays(rays, regions, length = sideLen)
        plot_k(np.arange(counter+1),kVals, ktitle)
        if ngroup == 10:
            energy_groups = [0.0,0.058,0.14,0.28,0.625,4.0,10.0,40.0,5530.0,821e3,20e6]
            plot_flux(energy_groups, regions)

    return k, regions


