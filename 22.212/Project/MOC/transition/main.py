NGROUP = 10
import numpy as np
import time
from numpy.random import random_sample as rand
from math import pi
from surface import *
from tools import *
from plotting import *
from region import *
from ray import *
from physics import *

np.random.seed(42)

def main(n_rays, surfaces, regions, sideLen, ngroup, plot=False, maxRayDist=100.0, deadzone=10):
    start = time.perf_counter()
    rays = []
    allActiveDists = 0.0

    rays = [drawRay(Ray(sideLen,maxRayDist), surfaces, regions, deadzone) for i in range(n_rays)]
    allActiveDists = sum([ray.active_length for ray in rays])

    for region in regions:
        region.vol = region.activeDist/allActiveDists
    
    counter = 0
    oldFissSrc, k = calc_q(regions, ngroup, 1, justStarting=True)
    kVals = [k]
    converged = False
    while not converged:
        normalize_phi(regions, ngroup)
        counter += 1
        print('Iterations: ', counter, ' k = ', k)

        updatePhi(rays, ngroup, regions)

        for reg in regions:
            SigT, vol = reg.mat.SigT, reg.vol
            reg.phi        = reg.tracks_phi/(vol*SigT*allActiveDists) + reg.q/SigT
            reg.tracks_phi = np.zeros([ngroup,])

        newFissSrc, k = calc_q(regions, ngroup, k, oldFissSrc)
        diffFissSrc = abs(newFissSrc-oldFissSrc)/oldFissSrc
        diffK       = abs((kVals[-1]-k)/kVals[-1])
        oldFissSrc = newFissSrc

        converged = diffK < 1e-5 and diffFissSrc < 1e-7
        kVals.append(k)
        
        if counter > 500: break
        
    print('k = ', k, ' after ', counter, 'iterations')
    end = time.perf_counter()
    elapsed_time = end - start

    segments = sum([len(ray.segments) for ray in rays])

    print('Elapsed time:               ', elapsed_time)
    print('Time per segment per group: ', elapsed_time/(segments*ngroup))

    if plot:
        ktitle ='k = '+str(k)+' Rays ='+str(n_rays)
        print('Plotting tracks')
        plot_from_rays(rays, regions, length = sideLen)
        plot_k(np.arange(counter+1),kVals, ktitle)
        if ngroup == 10:
            energy_groups = [0.0,0.058,0.14,0.28,0.625,4.0,10.0,40.0,5530.0,821e3,20e6]
            plot_flux(energy_groups, regions)

    return k, regions


