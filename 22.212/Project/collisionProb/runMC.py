import time
from numpy.random import random_sample as rand
from plotting import *
from ray import *
from physics import *
from materials import *

np.random.seed(42)

def runMC(n_rays, surfaces, regions, sideLen, ngroup, plot=False, 
            maxRayDist=100.0, deadzone=10, verbose=True, sphIter=0):

    start = time.perf_counter()

    ##########################################################################
    # DRAW NEUTRON 
    ##########################################################################

    pinHits = np.array([0]*9)
    modHits = 0
     
    circles = surfaces[0]
    if plot:
        ax = plt.gca()
        ax.cla()
        ax.set_xlim((0,sideLen))
        ax.set_ylim((0,sideLen))

        for c in circles:
            cPlot = plt.Circle((c.x0,c.y0),c.r,color='r')
            ax.add_artist(cPlot)
        colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
                  'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']




    numParticles = 1e4
    for i in range(numParticles):
        n = Neutron(sideLen,circles[0])
        if plot:
            ax.plot(n.r[0],n.r[1],'bo')

        counter = 0
        while True:
            n_next,distTraveled,regionID = advance(n,surfaces,regions)
            if plot:
                ax.plot(n_next.r[0],n_next.r[1],color=colors[counter%len(colors)],marker='o')
                ax.plot([n.r[0],n_next.r[0]],[n.r[1],n_next.r[1]],colors[counter%len(colors)]) 
            counter += 1
            #print(n.r,n.u)
            #print(n_next.r,n_next.u)
            #print()

            if weCollide(distTraveled,regionID,regions):
                #hits[regionID] += 1
                if regions[regionID].name == 'mod':
                    modHits += 1
                else:
                    pinHits[regionID] += 1
                break

            n.r = n_next.r
            n.u = n_next.u
    print(modHits/numParticles)
    print(pinHits/numParticles*100)
    assert(sum(pinHits)+modHits == numParticles)

    if plot:
        plt.show()
    






