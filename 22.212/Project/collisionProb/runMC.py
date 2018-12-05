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
    # DRAW ALL RAYS
    ##########################################################################
    #Ray(sideLen,maxRayDist), surfaces, regions, deadzone)


    circles = surfaces[0]
    ax = plt.gca()
    ax.cla()
    ax.set_xlim((0,sideLen))
    ax.set_ylim((0,sideLen))

    
    for c in circles:
        cPlot = plt.Circle((c.x0,c.y0),c.r,color='r')
        ax.add_artist(cPlot)


    n = Neutron(sideLen,circles[0])
    ax.plot([n.r[0],n.r[0]+n.u[0]],[n.r[1],n.r[1]+n.u[1]],'b') 
    ax.plot(n.r[0],n.r[1],'bo')

    n_next = advance(n,surfaces,regions)
    ax.plot(n_next.r[0],n_next.r[1],'go')


    plt.show()
    






