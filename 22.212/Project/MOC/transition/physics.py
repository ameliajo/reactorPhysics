import numpy as np
from math      import pi
from materials import *

def normPhiInAllRegs(regions, ngroup):
    phi_sum = sum([sum([reg.phi[g] for g in range(ngroup)]) for reg in regions])
    for region in regions:
        region.phi/= phi_sum# = [phi_g/phi_sum for phi_g in region.phi]



def updateQ(regions, ngroup, k, oldFissSrc=1,justStarting=False):
    newFissSrc = 0
    for reg in regions:
        mat = reg.mat
        regFissSrc = np.dot(mat.nuSigF,reg.phi)
        reg.q = np.matmul(mat.SigS,reg.phi) + mat.chi*regFissSrc/k
        newFissSrc += regFissSrc 

    k = newFissSrc/oldFissSrc 

    return (newFissSrc, 1) if justStarting else (newFissSrc/k,k)


def updatePhi(rays, ngroup, regions):
    for ray in rays:
        region = regions[ray.segments[0].region]
        mat = region.mat
        psi = region.q/(4*pi*mat.SigT)

        for segment in ray.segments:
            d = segment.d
            region = regions[segment.region]
            mat = region.mat
            q = region.q
            delta_psi = (psi - q/(4.0*pi*mat.SigT))*(1.0-np.exp(-mat.SigT*d))
            region.tracks_phi += 4.0*pi*delta_psi if segment.active else 0
            psi -= delta_psi
            
