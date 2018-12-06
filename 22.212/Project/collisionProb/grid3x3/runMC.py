import time
from neutron import *
from materials import *
import matplotlib.pyplot as plt

np.random.seed(42)

plot          = False
pitch         = 1.26
radius        = 0.39218
numParticles  = 1000



def runMC(pitch,radius,plot=False,numParticles=100):
    sideLen  = pitch*3

    regions,surfaces = makeGeometry(pitch,radius)

    start = time.perf_counter()

    pinHits, modHits = np.array([0]*9), 0
     
    circles = surfaces[0]
    if plot:
        ax = plt.gca(); ax.cla()
        ax.set_xlim((0,sideLen))
        ax.set_ylim((0,sideLen))

        for c in circles:
            cPlot = plt.Circle((c.x0,c.y0),c.r,color='r')
            ax.add_artist(cPlot)
        colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
                  'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']


    counter = 0
    for i in range(numParticles):
        n = Neutron(sideLen,circles[0])
        if plot: ax.plot(n.r[0],n.r[1],'bo')

        collided = False
        while not collided:
            if plot: oldr = n.r
            distTraveled,regionID = advance(n,surfaces,regions)
            if plot:
                ax.plot(n.r[0],n.r[1],color=colors[counter%len(colors)],marker='o')
                ax.plot([oldr[0],n.r[0]],[oldr[1],n.r[1]],colors[counter%len(colors)]) 
                counter += 1

            collided = weCollide(distTraveled,regions[regionID].mat)
            if collided and regions[regionID].name == 'mod':
                modHits += 1
            elif collided and regions[regionID].name == 'fuel':
                pinHits[regionID] += 1


    assert(sum(pinHits)+modHits == numParticles)

    if plot: plt.show()
    
    end = time.perf_counter()
    elapsed_time = end - start

    print('Elapsed time:', '%.4f'%elapsed_time)
    return pinHits/numParticles




collisionProb = runMC(pitch,radius,False,numParticles)
print(collisionProb*100)

