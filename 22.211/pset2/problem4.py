import numpy as np
import matplotlib.pyplot as plt
from math import log

def nslowing(nBins,nNeutrons):
    eBins    = np.linspace(0,999,nBins)
    lethBins = np.logspace(0,2.99,nBins)
    e_freq    = np.zeros(nBins-1)
    leth_freq = e_freq
    collision_dist = {1: e_freq,2: e_freq,3: e_freq}
    
    for i in range(int(nNeutrons)):
        energy = 1000 
        collisions = 0
        while energy > 1:
            e_freq_step    = np.histogram(energy,eBins)[0]
            leth_freq_step = np.histogram(energy,lethBins)[0]
            e_freq    += e_freq_step/20.0
            leth_freq += leth_freq_step/20.0
            if collisions in [1,2,3]:
                collision_dist[collisions] = collision_dist[collisions]+e_freq_step/20.0
            energy = np.random.uniform(0,1)*energy
            collisions += 1
    
    e_freq    = np.concatenate([[0],e_freq])   /nNeutrons
    leth_freq = np.concatenate([[0],leth_freq])/nNeutrons
    for i in range(3):
        collision_dist[i+1] = np.concatenate([[0],collision_dist[i+1]])/nNeutrons
    return eBins, lethBins, e_freq, leth_freq, collision_dist



if __name__ == '__main__':
    eBins, lethBins, e_freq, leth_freq, collision_dist = nslowing(nBins=100,nNeutrons=1e5)
    plt.title('Equal Energy Bins'); plt.ylabel('Normalized Flux (#/src neutron/cm^2)')
    plt.xlabel('Energy [eV]');     plt.step(eBins,e_freq); plt.show()

    plt.step(lethBins,leth_freq); 
    plt.title('Equal Lethargy Bins'); plt.ylabel('Normalized Flux (#/src neutron/cm^2)')
    plt.xlabel('Energy [eV]');        plt.yscale('log'); plt.show()
    
    for i in range(3):
        plt.step(lethBins,collision_dist[i+1])
    plt.legend(['1 Collison','2 Collisions','3 Collisions'])
    plt.title('Energy Distribution after collisions (Equal Lethargy Bins)')
    plt.xlabel('Energy [eV]'); plt.ylabel('Normalized Flux (#/src neutron/cm^2)')
    plt.yscale('log'); plt.show()


