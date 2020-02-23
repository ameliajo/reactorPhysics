import numpy as np
import matplotlib.pyplot as plt
from math import log
def nslowing(nBins):
    nNeutrons = int(1e4)
    eBins    = np.linspace(0,999,nBins)
    lethBins = np.logspace(0,2.99,nBins)
    e_freq    = np.zeros(nBins-1)
    leth_freq = e_freq
    collision_dist = {1: e_freq,2: e_freq,3: e_freq}
    
    for i in range(nNeutrons):
        energy = 1000 
        collisions = 0
        while energy > 1:
            e_freq_step    = np.histogram(energy,eBins)[0]
            leth_freq_step = np.histogram(energy,lethBins)[0]
            e_freq    += e_freq_step
            leth_freq += leth_freq_step
            if collisions in [1,2,3]:
                collision_dist[collisions] = collision_dist[collisions]+e_freq_step
            energy = np.random.uniform(0,1)*energy
            collisions += 1
    
    e_freq    = np.concatenate([[0],e_freq])   /nNeutrons
    leth_freq = np.concatenate([[0],leth_freq])/nNeutrons
    for i in range(3):
        collision_dist[i+1] = np.concatenate([[0],collision_dist[i+1]])/nNeutrons
    return eBins, lethBins, e_freq, leth_freq, collision_dist



if __name__ == '__main__':
    eBins, lethBins, e_freq, leth_freq, collision_dist = nslowing(nBins=100)
    ax1 = plt.subplot(3,1,1)
    plt.step(eBins,e_freq)
    plt.title('Equal Energy Bins')

    ax2 = plt.subplot(3,1,2)
    ax2.set_xscale("log", nonposx='clip')
    plt.step(lethBins,leth_freq)
    plt.ylabel('Normalized Flux (#/src neutron/cm^2)')
    plt.title('Equal Lethargy Bins')
    
    ax3=plt.subplot(3,1,3)
    ax3.set_xscale("log", nonposx='clip')
    ax3.set_yscale("log", nonposy='clip')
    for i in range(3):
        plt.step(lethBins,collision_dist[i+1])
    plt.legend(['1 Collison','2 Collisions','3 Collisions'])
    plt.title('Energy Distribution after collisions (Equal Lethargy)')
    plt.xlabel('Energy [eV]')
    plt.show()


